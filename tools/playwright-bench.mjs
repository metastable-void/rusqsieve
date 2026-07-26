// Real-browser benchmark for the published Web Worker architecture.
//
// Playwright stays an external development tool so the zero-dependency crate
// and browser artifact do not acquire a Node package manifest:
//
//   PLAYWRIGHT_MODULE=/tmp/rusqsieve-playwright/node_modules/playwright \
//   PLAYWRIGHT_BROWSERS_PATH=/tmp/rusqsieve-playwright-browsers \
//     node tools/playwright-bench.mjs URL N P Q [WORKERS]

import { createRequire } from "node:module";

const [url, decimal, expectedP, expectedQ, workersText = "8"] = process.argv.slice(2);
if (!url || !decimal || !expectedP || !expectedQ) {
  throw new Error("usage: playwright-bench.mjs URL N P Q [WORKERS]");
}
const workers = Number(workersText);
if (!Number.isInteger(workers) || workers < 1 || workers > 48) {
  throw new Error("WORKERS must be an integer from 1 through 48");
}

const require = createRequire(import.meta.url);
const playwrightModule = process.env.PLAYWRIGHT_MODULE || "playwright";
const { chromium } = require(playwrightModule);
const browser = await chromium.launch({ headless: true });
const context = await browser.newContext();
await context.addInitScript((count) => {
  Object.defineProperty(Navigator.prototype, "hardwareConcurrency", {
    configurable: true,
    get: () => count,
  });
}, workers);
const page = await context.newPage();
const browserErrors = [];
page.on("console", (message) => {
  if (message.type() !== "error") return;
  browserErrors.push(message.text());
  console.error(`browser console: ${message.text()}`);
});
page.on("pageerror", (error) => {
  browserErrors.push(error.message);
  console.error(`browser pageerror: ${error.message}`);
});

try {
  await page.goto(url, { waitUntil: "networkidle" });
  await page.locator("#go").waitFor({ state: "visible" });
  await page.waitForFunction(() => !document.querySelector("#go").disabled);
  const flavor = (await page.locator("#workers").textContent()).trim();
  if (!flavor.startsWith(`${workers} worker`)) {
    throw new Error(`worker override failed: ${flavor}`);
  }

  await page.locator("#input").evaluate((input, value) => {
    input.value = value;
    input.dispatchEvent(new Event("input", { bubbles: true }));
  }, decimal);
  const entered = await page.locator("#input").inputValue();
  if (entered !== decimal) {
    throw new Error(`input fill mismatch: got ${JSON.stringify(entered)}`);
  }
  await page.evaluate(() => {
    window.__rusqsieveBenchStart = performance.now();
    window.__rusqsieveBenchEvents = [];
    const status = document.querySelector("#status");
    new MutationObserver(() => {
      window.__rusqsieveBenchEvents.push({
        milliseconds: performance.now() - window.__rusqsieveBenchStart,
        text: status.textContent,
      });
    }).observe(status, { childList: true, subtree: true, characterData: true });
  });
  const runStarted = performance.now();
  await page.locator("#go").click();
  const progress = setInterval(async () => {
    try {
      const status = (await page.locator("#status").textContent()).trim();
      console.error(
        `playwright progress ${((performance.now() - runStarted) / 1000).toFixed(1)}s: ${status}`,
      );
    } catch {
      // The final navigation/close path may race this diagnostic poll.
    }
  }, 5000);
  try {
    await page.waitForFunction(
      () => {
        const status = document.querySelector("#status").textContent;
        return (
          status === "Done." ||
          status.startsWith("Error:") ||
          status === "Enter a positive whole number."
        );
      },
      undefined,
      { timeout: 30 * 60 * 1000 },
    );
  } finally {
    clearInterval(progress);
  }
  const finalStatus = await page.locator("#status").textContent();
  if (finalStatus !== "Done.") {
    throw new Error(`frontend failed: ${finalStatus}`);
  }

  const result = await page.evaluate(() => ({
    elapsedMilliseconds: performance.now() - window.__rusqsieveBenchStart,
    events: window.__rusqsieveBenchEvents,
    meta: document.querySelector("#result .meta")?.textContent || "",
    factors: document.querySelector("#result code")?.textContent || "",
  }));
  const expected = [BigInt(expectedP), BigInt(expectedQ)]
    .sort((a, b) => (a < b ? -1 : a > b ? 1 : 0))
    .join(" * ");
  if (!result.meta.includes("✓ verified") || result.factors !== expected) {
    throw new Error(
      `factor verification mismatch: meta=${result.meta} factors=${result.factors} expected=${expected}`,
    );
  }

  const firstSieve = result.events.find((event) => event.text.startsWith("Quadratic sieve:"));
  const linearAlgebra = result.events.find((event) => event.text.startsWith("Linear algebra"));
  const total = result.elapsedMilliseconds / 1000;
  const firstRelation = firstSieve ? firstSieve.milliseconds / 1000 : Number.NaN;
  const sieveEnd = linearAlgebra ? linearAlgebra.milliseconds / 1000 : Number.NaN;
  const linearAlgebraSeconds = linearAlgebra
    ? (result.elapsedMilliseconds - linearAlgebra.milliseconds) / 1000
    : Number.NaN;
  console.log(
    [
      `playwright-chromium ${flavor}`,
      `bits=${BigInt(decimal).toString(2).length}`,
      `total=${total.toFixed(3)}s`,
      `first_relation=${firstRelation.toFixed(3)}s`,
      `sieve_end=${sieveEnd.toFixed(3)}s`,
      `linalg_extract=${linearAlgebraSeconds.toFixed(3)}s`,
    ].join(" "),
  );
  console.log(result.meta);
  if (browserErrors.length) {
    throw new Error(`browser console errors: ${browserErrors.join(" | ")}`);
  }
} finally {
  await browser.close();
}
