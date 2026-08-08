//! Raw, versioned C ABI for `wasm32-unknown-unknown`.
use crate::engine::{self, EngineContext, EngineJob, EngineSession};
use crate::{Natural, PARTS};
use std::alloc::{Layout, alloc, dealloc};
use std::cell::RefCell;

// 3 adds `qs_max_siqs_bits` and `qs_coord_family_budget`. The glue does not merely tolerate these
// — it reads the sieve range and the family budget from them instead of hard-coding either — so a
// v2 module paired with current glue would fault on a missing export rather than degrade. The
// version guard has to reject that pairing.
//
// 4 adds `qs_rho`. The frontend's rho workers depend on it rather than falling back to their
// BigInt implementation when it is absent, for the same reason: a silent eightfold slowdown in the
// one stage that decides whether a wide multi-factor composite finishes is not a degradation worth
// shipping quietly.
//
// 5 adds `qs_ecm` and its two default-bound queries. The frontend runs curves on any composite the
// sieve refuses, so a module without them would silently lose the only stage that can factor a wide
// number with a 25-digit factor.
const ABI_VERSION: u32 = 5;
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
thread_local! {static BUFFERS:RefCell<Registry<Box<[u8]>>>=const{RefCell::new(Registry::new())};
    static COORDS: RefCell<Registry<EngineSession>> = const { RefCell::new(Registry::new()) };
    static WORKERS: RefCell<Registry<EngineContext>> = const { RefCell::new(Registry::new()) };
}
fn memory_bytes() -> usize {
    core::arch::wasm32::memory_size(0) * 65536
}
fn input(pointer: u32, length: u32) -> Option<Vec<u8>> {
    let p = pointer as usize;
    let n = length as usize;
    if n == 0 || n > MAX_PACKET || p.checked_add(n)? > memory_bytes() {
        return None;
    }
    // SAFETY: the explicit linear-memory bound check proves that this range is
    // readable. Copying immediately prevents later registry mutations from
    // invalidating or aliasing a borrowed caller buffer.
    Some(unsafe { core::slice::from_raw_parts(p as *const u8, n) }.to_vec())
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
/// Widest composite the sieve accepts. The coordinator refuses anything above this, so the glue
/// can reject it up front with a specific message instead of reporting a generic setup failure.
#[unsafe(no_mangle)]
pub extern "C" fn qs_max_siqs_bits() -> u32 {
    engine::MAX_SIQS_BITS as u32
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
pub unsafe extern "C" fn qs_dealloc(pointer: u32, size: u32, align: u32) {
    let Ok(layout) = Layout::from_size_align(size as usize, align as usize) else {
        return;
    };
    if size != 0 && pointer != 0 {
        // SAFETY: the caller must pass the same pointer and layout previously
        // returned by `qs_alloc`, and must not have freed it already.
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
    let text = core::str::from_utf8(&bytes).ok()?;
    WasmNatural::from_decimal(text).ok()
}

/// Prepare a deterministic worker sieve context for the composite `n`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn qs_worker_prepare(n_pointer: u32, n_length: u32) -> u32 {
    let Some(n) = parse_decimal(n_pointer, n_length) else {
        return 0;
    };
    let Ok(ctx) = engine::prepare(n, &crate::factor::FactorTuning::default()) else {
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
pub unsafe extern "C" fn qs_coord_new(n_pointer: u32, n_length: u32) -> u32 {
    let Some(n) = parse_decimal(n_pointer, n_length) else {
        return 0;
    };
    let Ok(ctx) = engine::prepare(n, &crate::factor::FactorTuning::default()) else {
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
/// Polynomial families this session will issue before the relation budget is spent.
///
/// The budget scales with input width, so a scheduler that assigns family numbers itself — as the
/// browser coordinator does — must read it here rather than hard-coding a constant.
#[unsafe(no_mangle)]
pub extern "C" fn qs_coord_family_budget(session: u32) -> u32 {
    COORDS.with(|r| {
        r.borrow()
            .get(session)
            .map_or(0, |s| u32::try_from(s.family_budget()).unwrap_or(u32::MAX))
    })
}
/// Ingest a worker's `qs_worker_sieve` payload; returns the new relation count.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn qs_coord_submit(session: u32, pointer: u32, length: u32) -> u32 {
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

// ---------------------------------------------------------------------------
// Deep Pollard-Brent for the browser's rho workers.
//
// The frontend peels cheap factors on the main thread with BigInt, which is fine for an
// opportunistic peel and hopeless for a real search: measured under Node, BigInt runs this loop at
// about an eighth of the rate the same algorithm reaches here over Montgomery-encoded limbs. A
// composite the sieve cannot help with — one above `qs_max_siqs_bits`, or a cofactor an earlier
// split already proved unbalanced — needs tens of millions of iterations, and that difference is
// what decides whether a 48-bit factor is found or the number is handed to a sieve that will not
// finish.
//
// Cancellation is the worker's `terminate()`, so nothing here polls: the call runs to its budget
// and returns. That is also why the budget is a parameter rather than a policy baked in here — the
// glue sizes it per width, exactly as `engine::rho_budget` does natively.
// ---------------------------------------------------------------------------

/// Search for a factor of the decimal composite `n` with the elliptic curve method.
///
/// ECM is the tool for a medium-size factor — 20 to 30 digits — which is the shape Pollard-Brent
/// cannot reach and the sieve either pays for by the size of `n` or refuses outright. The engine
/// runs it by itself on a composite the sieve refuses or one already known to be unbalanced; this
/// export exists so the browser can do the same from a worker. Its cost
/// depends on the size of the factor rather than of the input, so a wide composite with a 25-digit
/// factor is ordinary work here and impossible for the other two.
///
/// `b1`/`b2` are the stage bounds and `curves` the number of curves to try; `seed` selects the
/// `σ` sequence, so the same arguments reproduce the same run. Zero bounds take the defaults for
/// the composite's width.
///
/// Returns a packet (kind 13) whose payload is the factor as `PARTS * 8` little-endian bytes, or a
/// packet with an empty payload when every curve was exhausted. Returns 0 only for a modulus that
/// cannot be parsed or is below 3. Cancellation is the worker's `terminate()`, as for `qs_rho`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn qs_ecm(
    n_pointer: u32,
    n_length: u32,
    b1: u32,
    b2: u32,
    curves: u32,
    seed: u32,
) -> u32 {
    let Some(n) = parse_decimal(n_pointer, n_length) else {
        return 0;
    };
    if n < WasmNatural::from_u64(3) {
        return 0;
    }
    let defaults = crate::ecm::EcmParams::for_composite(n.bit_len());
    let params = crate::ecm::EcmParams {
        b1: if b1 == 0 { defaults.b1 } else { u64::from(b1) },
        b2: if b2 == 0 { defaults.b2 } else { u64::from(b2) },
        curves: if curves == 0 { defaults.curves } else { curves },
    };
    let found = crate::ecm::factor(&n, params, u64::from(seed), || true);
    let mut payload = Vec::new();
    if let Ok(Some(factor)) = found {
        payload.reserve(PARTS * 8);
        for limb in factor.as_parts() {
            payload.extend_from_slice(&limb.to_le_bytes());
        }
    }
    packet(13, &payload)
}

/// Report the default ECM bounds for a composite of `bits` bits, so the glue can size a run and
/// report it without duplicating the schedule.
#[unsafe(no_mangle)]
pub extern "C" fn qs_ecm_default_b1(bits: u32) -> u32 {
    u32::try_from(crate::ecm::EcmParams::for_composite(bits as usize).b1).unwrap_or(u32::MAX)
}

/// Companion to [`qs_ecm_default_b1`].
#[unsafe(no_mangle)]
pub extern "C" fn qs_ecm_default_curves(bits: u32) -> u32 {
    crate::ecm::EcmParams::for_composite(bits as usize).curves
}

/// Search for a factor of the decimal composite `n` with a bounded Pollard-Brent.
///
/// `first_constant` and `constant_count` select the polynomial constants `y^2 + c` this call walks,
/// so a pool of workers given disjoint ranges runs that many independent walks over the same
/// modulus and the first collision wins. `budget` is the total iterations across those walks.
///
/// Returns a packet (kind 12) whose payload is the factor as `PARTS * 8` little-endian bytes, or a
/// packet with an empty payload when the budget was spent without a split. Returns 0 only for an
/// unusable request: an unparseable modulus, or one below 3.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn qs_rho(
    n_pointer: u32,
    n_length: u32,
    budget: u32,
    first_constant: u32,
    constant_count: u32,
) -> u32 {
    let Some(n) = parse_decimal(n_pointer, n_length) else {
        return 0;
    };
    // Below 3 there is nothing to split, and a modulus of 0 or 1 has no Montgomery form.
    if n < WasmNatural::from_u64(3) {
        return 0;
    }
    let first = (first_constant.max(1)) as u64;
    // The cap keeps `first + count` from wrapping and matches what a pool could plausibly use;
    // beyond a handful of constants the budget per walk is what limits reach anyway.
    let count = constant_count.clamp(1, 64) as u64;
    let last = first.saturating_add(count - 1);
    let found = engine::pollard_brent_natural(&n, (budget as u64).max(1), first..=last, || true);
    let mut payload = Vec::new();
    if let Ok(Some(factor)) = found {
        payload.reserve(PARTS * 8);
        for limb in factor.as_parts() {
            payload.extend_from_slice(&limb.to_le_bytes());
        }
    }
    packet(12, &payload)
}
