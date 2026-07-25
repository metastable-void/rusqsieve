// Runs web/coordinator.js on a node worker thread, giving it the `self.onmessage` /
// `self.postMessage` surface a browser Worker provides. Mirrors node-worker-shim.mjs,
// which does the same for the sieve worker; together they let the browser's
// coordinator-in-its-own-Worker architecture be exercised without a browser.
import { parentPort } from "node:worker_threads";

const pending = [];
globalThis.self = {
  onmessage: null,
  postMessage(message, transfer) {
    parentPort.postMessage(message, transfer);
  },
};

parentPort.on("message", (data) => {
  if (self.onmessage) self.onmessage({ data });
  else pending.push(data);
});

await import("../web/coordinator.js");
for (const data of pending.splice(0)) self.onmessage({ data });
