//! Raw, versioned C ABI for `wasm32-unknown-unknown`.
use crate::engine::{self, EngineContext, EngineJob, EngineSession};
use crate::{FactorConfig, FactorSession, LocalWorkBudget, Natural, PARTS};
use std::alloc::{Layout, alloc, dealloc};
use std::cell::RefCell;

const ABI_VERSION: u32 = 1;
const MAX_PACKET: usize = 16 * 1024 * 1024;
type WasmNatural = Natural;
struct Slot<T> {
    generation: u16,
    value: Option<T>,
}
struct Registry<T> {
    slots: Vec<Slot<T>>,
}
impl<T> Registry<T> {
    const fn new() -> Self {
        Self { slots: Vec::new() }
    }
    fn insert(&mut self, value: T) -> u32 {
        for (i, s) in self.slots.iter_mut().enumerate() {
            if s.value.is_none() {
                s.value = Some(value);
                return ((s.generation as u32) << 16) | (i as u32 + 1);
            }
        }
        let i = self.slots.len();
        if i >= u16::MAX as usize {
            return 0;
        }
        self.slots.push(Slot {
            generation: 1,
            value: Some(value),
        });
        ((1u32) << 16) | (i as u32 + 1)
    }
    fn get(&self, h: u32) -> Option<&T> {
        let i = (h & 0xffff).checked_sub(1)? as usize;
        let g = (h >> 16) as u16;
        let s = self.slots.get(i)?;
        if s.generation == g {
            s.value.as_ref()
        } else {
            None
        }
    }
    fn get_mut(&mut self, h: u32) -> Option<&mut T> {
        let i = (h & 0xffff).checked_sub(1)? as usize;
        let g = (h >> 16) as u16;
        let s = self.slots.get_mut(i)?;
        if s.generation == g {
            s.value.as_mut()
        } else {
            None
        }
    }
    fn remove(&mut self, h: u32) {
        let i = match (h & 0xffff).checked_sub(1) {
            Some(v) => v as usize,
            None => return,
        };
        let g = (h >> 16) as u16;
        if let Some(s) = self.slots.get_mut(i)
            && s.generation == g
        {
            s.value = None;
            s.generation = s.generation.wrapping_add(1).max(1)
        }
    }
}
thread_local! {static SESSIONS:RefCell<Registry<FactorSession<PARTS>>>=const{RefCell::new(Registry::new())};static BUFFERS:RefCell<Registry<Box<[u8]>>>=const{RefCell::new(Registry::new())};
    static COORDS: RefCell<Registry<EngineSession>> = const { RefCell::new(Registry::new()) };
    static WORKERS: RefCell<Registry<EngineContext>> = const { RefCell::new(Registry::new()) };
}
fn memory_bytes() -> usize {
    core::arch::wasm32::memory_size(0) * 65536
}
fn input(pointer: u32, length: u32) -> Option<&'static [u8]> {
    let p = pointer as usize;
    let n = length as usize;
    if n == 0 || n > MAX_PACKET || p.checked_add(n)? > memory_bytes() {
        return None;
    }
    Some(unsafe { core::slice::from_raw_parts(p as *const u8, n) })
}
fn packet(kind: u16, payload: &[u8]) -> u32 {
    if payload.len() > MAX_PACKET {
        return 0;
    }
    let mut v = Vec::with_capacity(12 + payload.len());
    v.extend_from_slice(b"QSV1");
    v.extend_from_slice(&kind.to_le_bytes());
    v.extend_from_slice(&1u16.to_le_bytes());
    v.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    v.extend_from_slice(payload);
    BUFFERS.with(|r| r.borrow_mut().insert(v.into_boxed_slice()))
}

#[unsafe(no_mangle)]
pub extern "C" fn qs_abi_version() -> u32 {
    ABI_VERSION
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_alloc(size: u32, align: u32) -> u32 {
    let Ok(layout) = Layout::from_size_align(size as usize, align as usize) else {
        return 0;
    };
    if size == 0 {
        return align;
    }
    let p = unsafe { alloc(layout) };
    p as u32
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_dealloc(pointer: u32, size: u32, align: u32) {
    let Ok(layout) = Layout::from_size_align(size as usize, align as usize) else {
        return;
    };
    if size != 0 && pointer != 0 {
        unsafe { dealloc(pointer as *mut u8, layout) }
    }
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_buffer_pointer(handle: u32) -> u32 {
    BUFFERS.with(|r| r.borrow().get(handle).map_or(0, |b| b.as_ptr() as u32))
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_buffer_length(handle: u32) -> u32 {
    BUFFERS.with(|r| {
        r.borrow()
            .get(handle)
            .and_then(|b| u32::try_from(b.len()).ok())
            .unwrap_or(0)
    })
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_buffer_free(handle: u32) {
    BUFFERS.with(|r| r.borrow_mut().remove(handle))
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_new(
    input_pointer: u32,
    input_length: u32,
    _config_pointer: u32,
    _config_length: u32,
) -> u32 {
    let Some(bytes) = input(input_pointer, input_length) else {
        return 0;
    };
    let Ok(text) = core::str::from_utf8(bytes) else {
        return 0;
    };
    let Ok(n) = WasmNatural::from_decimal(text) else {
        return 0;
    };
    let Ok(s) = FactorSession::new(n, FactorConfig::default()) else {
        return 0;
    };
    SESSIONS.with(|r| r.borrow_mut().insert(s))
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_free(session: u32) {
    SESSIONS.with(|r| r.borrow_mut().remove(session))
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_phase(session: u32) -> u32 {
    SESSIONS.with(|r| {
        r.borrow()
            .get(session)
            .map_or(u32::MAX, |s| s.phase() as u32)
    })
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_advance_local(session: u32, _p: u32, _n: u32) -> i32 {
    SESSIONS.with(|r| match r.borrow_mut().get_mut(session) {
        Some(s) => match s.advance_local(LocalWorkBudget::default()) {
            Ok(crate::AdvanceOutcome::Complete) => 1,
            Ok(_) => 0,
            Err(_) => -1,
        },
        None => -2,
    })
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_export_context(_session: u32) -> u32 {
    packet(2, &[])
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_take_jobs(_session: u32, _maximum_jobs: u32) -> u32 {
    packet(3, &[])
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_submit(_session: u32, _pointer: u32, _length: u32) -> i32 {
    -1
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_take_factors(session: u32) -> u32 {
    SESSIONS.with(|r| {
        let mut r = r.borrow_mut();
        let Some(i) = ((session & 0xffff).checked_sub(1)).map(|x| x as usize) else {
            return 0;
        };
        let Some(slot) = r.slots.get_mut(i) else {
            return 0;
        };
        if slot.generation != (session >> 16) as u16 {
            return 0;
        }
        let Some(s) = slot.value.take() else { return 0 };
        match s.take_factors() {
            Ok(fs) => {
                let mut text = String::new();
                for (p, e) in fs.iter() {
                    text.push_str(&p.to_string());
                    text.push(':');
                    text.push_str(&e.get().to_string());
                    text.push('\n')
                }
                packet(4, text.as_bytes())
            }
            Err(_) => 0,
        }
    })
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_error(_session: u32) -> u32 {
    packet(5, &[])
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_session_progress(session: u32) -> u32 {
    SESSIONS.with(|r| {
        r.borrow().get(session).map_or(0, |s| {
            let p = s.progress();
            let mut v = Vec::new();
            v.extend_from_slice(&p.revision.to_le_bytes());
            v.extend_from_slice(&(p.phase as u32).to_le_bytes());
            v.extend_from_slice(&p.amount.completed.to_le_bytes());
            packet(6, &v)
        })
    })
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_context_import(_pointer: u32, _length: u32) -> u32 {
    0
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_context_free(_context: u32) {}
#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_execute(_context: u32, _pointer: u32, _length: u32) -> u32 {
    0
}

// ---------------------------------------------------------------------------
// Parallel SIQS protocol (engine-based) for the Web-Worker demo.
//
// A worker rebuilds the *deterministic* sieve context with `qs_worker_prepare`
// (same input → same factor base, so no context needs to be serialized) and
// sieves a stripe of polynomial families with `qs_worker_sieve`. The coordinator
// (`qs_coord_*`) accumulates the serialized relations and runs the linear algebra.
// ---------------------------------------------------------------------------

fn parse_decimal(pointer: u32, length: u32) -> Option<WasmNatural> {
    let bytes = input(pointer, length)?;
    let text = core::str::from_utf8(bytes).ok()?;
    WasmNatural::from_decimal(text).ok()
}

/// Prepare a deterministic worker sieve context for the composite `n`.
#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_prepare(n_pointer: u32, n_length: u32) -> u32 {
    let Some(n) = parse_decimal(n_pointer, n_length) else {
        return 0;
    };
    let Ok(ctx) = engine::prepare(n) else {
        return 0;
    };
    WORKERS.with(|r| r.borrow_mut().insert(ctx))
}
/// Sieve polynomial families `[family_first, family_first + count)`; returns a buffer
/// handle to `count` concatenated serialized family results (`[count:u32][len:u32,bytes]…`).
#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_sieve(context: u32, family_first: u32, count: u32) -> u32 {
    WORKERS.with(|r| {
        let reg = r.borrow();
        let Some(ctx) = reg.get(context) else {
            return 0;
        };
        let count = count.min(4096);
        let mut payload = Vec::new();
        payload.extend_from_slice(&count.to_le_bytes());
        for k in 0..count {
            let job = EngineJob {
                family: (family_first + k) as u64,
            };
            let bytes = engine::execute(ctx, job).to_bytes();
            payload.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
            payload.extend_from_slice(&bytes);
        }
        packet(10, &payload)
    })
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_worker_free(context: u32) {
    WORKERS.with(|r| r.borrow_mut().remove(context))
}

/// Create a coordinator collecting relations for the composite `n`.
#[unsafe(no_mangle)]
pub extern "C" fn qs_coord_new(n_pointer: u32, n_length: u32) -> u32 {
    let Some(n) = parse_decimal(n_pointer, n_length) else {
        return 0;
    };
    let Ok(ctx) = engine::prepare(n) else {
        return 0;
    };
    COORDS.with(|r| r.borrow_mut().insert(EngineSession::new(ctx)))
}
/// Relation target needed before the coordinator can extract a factor.
#[unsafe(no_mangle)]
pub extern "C" fn qs_coord_target(session: u32) -> u32 {
    COORDS.with(|r| r.borrow().get(session).map_or(0, |s| s.target() as u32))
}
/// Relations collected so far.
#[unsafe(no_mangle)]
pub extern "C" fn qs_coord_relations(session: u32) -> u32 {
    COORDS.with(|r| r.borrow().get(session).map_or(0, |s| s.relations() as u32))
}
/// Ingest a worker's `qs_worker_sieve` payload; returns the new relation count.
#[unsafe(no_mangle)]
pub extern "C" fn qs_coord_submit(session: u32, pointer: u32, length: u32) -> u32 {
    let Some(bytes) = input(pointer, length) else {
        return 0;
    };
    COORDS.with(|r| {
        let mut reg = r.borrow_mut();
        let Some(s) = reg.get_mut(session) else {
            return 0;
        };
        if bytes.len() >= 4 {
            let count = u32::from_le_bytes(bytes[0..4].try_into().unwrap());
            let mut o = 4usize;
            for _ in 0..count {
                if o + 4 > bytes.len() {
                    break;
                }
                let len = u32::from_le_bytes(bytes[o..o + 4].try_into().unwrap()) as usize;
                o += 4;
                if o + len > bytes.len() {
                    break;
                }
                s.submit_bytes(&bytes[o..o + len]);
                o += len;
            }
        }
        s.relations() as u32
    })
}
/// Run the linear algebra and return a nontrivial factor as a `PARTS * 8`-byte
/// little-endian `Natural` payload, or 0 if extraction failed (needs more relations).
#[unsafe(no_mangle)]
pub extern "C" fn qs_coord_extract(session: u32) -> u32 {
    COORDS.with(|r| {
        let reg = r.borrow();
        let Some(s) = reg.get(session) else {
            return 0;
        };
        match s.extract_factor() {
            Ok(d) => {
                let mut payload = Vec::with_capacity(PARTS * 8);
                for limb in d.as_parts() {
                    payload.extend_from_slice(&limb.to_le_bytes());
                }
                packet(11, &payload)
            }
            Err(_) => 0,
        }
    })
}
#[unsafe(no_mangle)]
pub extern "C" fn qs_coord_free(session: u32) {
    COORDS.with(|r| r.borrow_mut().remove(session))
}
