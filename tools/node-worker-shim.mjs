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

await import("../web/worker.js");
for (const data of pending.splice(0)) self.onmessage({ data });
