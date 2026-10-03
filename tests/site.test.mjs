import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const app = fs.readFileSync(path.join(root, 'js/app.js'), 'utf8');

test('local scripts exist and stay relative for the /iron_oxydation/ subpath', () => {
  const srcs = [...html.matchAll(/<script src="([^"]+)"/g)].map((m) => m[1]);
  assert.deepEqual(srcs, ['js/vendor/mathjax/tex-svg.js', 'js/app.js']);
  for (const src of srcs) assert.ok(fs.existsSync(path.join(root, src)), src);
  assert.equal(/(?:src|href)="\//.test(html), false, 'no root-absolute paths');
});

test('equations are TeX source for MathJax, not a stale rendered copy', () => {
  assert.equal(/<mjx-container|MJX-SVG/.test(html), false);
  const display = (html.match(/\\\[/g) || []).length;
  const inline = (html.match(/\\\(/g) || []).length;
  // 71 equations plus the delimiter declarations in the MathJax config.
  assert.equal(display + inline, 73);
});

test('no claude.ai frame runtime and no network calls are left in the page', () => {
  for (const text of [html, app]) {
    assert.equal(/__FRAME|frame-runtime|claudeusercontent/.test(text), false);
    assert.equal(/\bfetch\(|XMLHttpRequest|window\.claude/.test(text), false);
  }
});

test('app script parses, JS-built rows start empty, and every element it looks up is in the page', () => {
  new vm.Script(app);
  for (const id of ['forms', 'metals', 'presets']) assert.ok(html.includes(`id="${id}"></div>`), id);
  for (const [, id] of app.matchAll(/\$\('([a-zA-Z]+)'\)/g)) assert.ok(html.includes(`id="${id}"`), id);
});
