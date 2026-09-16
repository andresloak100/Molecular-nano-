"use strict";

// Local evidence is data only. Imported names and text never become markup.
const $ = (id) => document.getElementById(id);
const state = { catalog: null, structure: null, tab: "structure", campaign: null, pose: null, endpoint: "initial", selected: [], loadToken: 0 };
const finite = (value) => typeof value === "number" && Number.isFinite(value);
const fmt = (value, digits = 3) => finite(value) ? value.toLocaleString("en-US", { minimumFractionDigits: digits, maximumFractionDigits: digits }) : "Not recorded";
const signed = (value, digits = 3) => finite(value) ? `${value > 0 ? "+" : ""}${fmt(value, digits)}` : "Not recorded";
const human = (value) => typeof value === "string" ? value.replaceAll("_", " ") : "Not recorded";
const yesNo = (value) => value === true ? "Yes" : value === false ? "No" : "Not recorded";
const array = (value) => Array.isArray(value) ? value : [];
const object = (value) => value && typeof value === "object" && !Array.isArray(value) ? value : {};
function node(tag, className, text) {
  const result = document.createElement(tag);
  if (className) result.className = className;
  if (text !== undefined) result.textContent = String(text);
  return result;
}
function append(parent, ...children) { parent.append(...children.filter(Boolean)); return parent; }
function statusPill(status) {
  const value = typeof status === "string" ? status : "not recorded";
  const known = ["completed", "pending", "failed", "error", "running"].includes(value) ? value : ["invalid evidence", "unreadable"].includes(value) ? "error" : "";
  return node("span", `status-pill ${known}`, human(value));
}
function sourceLink(record, label = "View source record ↗") {
  const href = record && record.source_href;
  if (typeof href !== "string" || !href.startsWith("/api/evidence?") || href.startsWith("//")) return null;
  const a = node("a", "text-link", label); a.href = href; a.target = "_blank"; a.rel = "noopener"; return a;
}
function detailList(pairs) {
  const dl = node("dl", "detail-list");
  pairs.forEach(([key, value]) => append(dl, append(node("div"), node("dt", "", key), node("dd", "", value ?? "Not recorded"))));
  return dl;
}
function metric(label, value) { return append(node("div"), node("span", "", label), node("strong", "", value)); }
function empty(title, body, icon = "◇") {
  return append(node("div", "content-card empty-state"), node("span", "empty-symbol", icon), node("h2", "", title), node("p", "", body));
}
function setError(message) { $("error-banner").hidden = !message; $("error-banner").textContent = message || ""; }
async function fetchJSON(path) {
  const response = await fetch(path, { headers: { Accept: "application/json" } });
  if (!response.ok) {
    let reason = `Local record could not be read (${response.status}).`;
    try { const body = await response.json(); if (typeof body.error === "string") reason = body.error; } catch (_) { /* Preserve HTTP context. */ }
    throw new Error(reason);
  }
  return response.json();
}
function distance(a, b, structure = state.structure) {
  if (!Number.isInteger(a) || !Number.isInteger(b)) return null;
  const p = structure?.positions?.[a], q = structure?.positions?.[b];
  return p && q ? Math.hypot(...p.map((v, i) => v - q[i])) : null;
}
function transfer(structure = state.structure) { return object(structure?.hydrogen_transfer || structure?.metadata?.hydrogen_transfer); }
function formula(symbols) {
  const counts = new Map(); symbols.forEach(s => counts.set(s, (counts.get(s) || 0) + 1));
  const keys = [...counts.keys()].sort((a, b) => (a === "C" ? -2 : a === "H" ? -1 : 0) - (b === "C" ? -2 : b === "H" ? -1 : 0) || a.localeCompare(b));
  return keys.map(s => `${s}${counts.get(s) > 1 ? counts.get(s) : ""}`).join("");
}
function structureMeta(structure) { return { ...object(structure?.recorded_result?.design?.metadata), ...object(structure?.metadata) }; }
function atomRole(index) {
  if (!state.structure) return "";
  const t = transfer(), meta = structureMeta(state.structure), labels = [];
  if (index === t.donor) labels.push("donor");
  if (index === t.hydrogen) labels.push("transferred H");
  if (index === t.acceptor) labels.push("acceptor");
  if (array(state.structure.fixed_indices).includes(index)) labels.push("fixed anchor");
  if (array(meta.tool_indices).includes(index)) labels.push("tool");
  else if (array(meta.substrate_indices).includes(index)) labels.push("substrate");
  return labels.join(" · ");
}
function updateMeasurements() {
  const [a, b] = state.selected;
  $("atom-a").value = a === undefined ? "" : String(a);
  $("atom-b").value = b === undefined ? "" : String(b);
  const value = distance(a, b);
  $("pair-distance").replaceChildren(document.createTextNode(finite(value) ? fmt(value, 4) + " " : "— "), node("span", "", "Å"));
  if (a === undefined) $("selection-detail").textContent = "Choose two atoms in the viewer or above.";
  else if (b === undefined) $("selection-detail").textContent = `${state.structure.symbols[a]} ${a} · ${atomRole(a) || "selected"}. Choose a second atom.`;
  else $("selection-detail").textContent = `${state.structure.symbols[a]} ${a} ↔ ${state.structure.symbols[b]} ${b} · Euclidean distance from recorded coordinates.`;
  viewer.draw();
}
function selectAtom(index) {
  if (state.selected.includes(index)) state.selected = state.selected.filter(i => i !== index);
  else if (state.selected.length < 2) state.selected.push(index);
  else state.selected = [index];
  updateMeasurements();
}
function updateStructureDetails() {
  const structure = state.structure, meta = structureMeta(structure), t = transfer();
  const symbols = array(structure.symbols), fixed = array(structure.fixed_indices);
  $("formula").textContent = formula(symbols) || "Not recorded";
  $("atom-count").textContent = symbols.length;
  $("anchor-count").textContent = fixed.length;
  const da = distance(t.donor, t.acceptor);
  $("apex-distance").textContent = finite(da) ? `${fmt(da, 3)} Å` : "Not recorded";
  $("structure-label").textContent = structure.label || "Recorded structure";
  $("structure-select").value = structure.id;
  $("bond-note").textContent = structure.bond_source || "Distance-inferred bonds are visual guides, not electronic bond orders.";
  const source = sourceLink(structure, "Source ↗");
  $("structure-source").hidden = !source;
  if (source) $("structure-source").href = source.href;
  const groups = [array(meta.substrate_indices).length, array(meta.tool_indices).length];
  const pairs = [["Geometry", human(meta.geometry_status || "not recorded")], ["Substrate / tool", groups.some(Boolean) ? `${groups[0]} / ${groups[1]} atoms` : "Not recorded"], ["Fixed atom indices", fixed.length ? fixed.join(", ") : "None recorded"], ["Charge / multiplicity", `${meta.total_charge ?? meta.charge ?? "—"} / ${meta.spin_multiplicity ?? meta.multiplicity ?? "—"}`], ["Boundary", meta.periodic === false ? "Nonperiodic" : meta.periodic === true ? "Periodic" : "Not recorded"]];
  $("geometry-details").replaceChildren(...detailList(pairs).children);
  $("geometry-scope").textContent = meta.description || structure.recorded_result?.design?.scope || "Only the saved geometry is shown. Atom positions are not changed by rotating the view.";
  for (const id of ["atom-a", "atom-b"]) {
    const select = $(id); select.replaceChildren(new Option("Select atom", ""));
    symbols.forEach((s, i) => select.add(new Option(`${s} ${i}`, String(i))));
  }
  const panel = $("transfer-details"); panel.replaceChildren();
  if ([t.donor, t.hydrogen, t.acceptor].every(i => Number.isInteger(i) && symbols[i])) {
    const chain = node("div", "transfer-chain");
    [[t.donor, "Donor"], [t.hydrogen, "Hydrogen"], [t.acceptor, "Acceptor"]].forEach(([index, name], i) => {
      if (i) chain.append(node("span", "transfer-arrow", "···"));
      const circle = node("button", "transfer-atom", `${symbols[index]} ${index}`);
      circle.title = `Select ${name.toLowerCase()} atom ${index}`; circle.addEventListener("click", () => selectAtom(index));
      chain.append(append(node("div"), circle, node("div", "transfer-caption", name)));
    });
    const distances = node("div", "transfer-distances");
    append(distances, append(node("div", "", "Donor–H"), node("strong", "", `${fmt(distance(t.donor, t.hydrogen), 3)} Å`)), append(node("div", "", "Acceptor–H"), node("strong", "", `${fmt(distance(t.acceptor, t.hydrogen), 3)} Å`)));
    append(panel, chain, distances);
  } else panel.append(node("p", "muted small", "No hydrogen-transfer assignment recorded for this structure."));
  const groupSelect = $("group-select");
  [...groupSelect.options].forEach(option => { option.disabled = option.value !== "all" && !array(meta[`${option.value}_indices`]).length; });
  if (groupSelect.selectedOptions[0]?.disabled) groupSelect.value = "all";
  $("coordinate-status").textContent = `${symbols.length} atoms · Actual coordinates · Å`;
  updateMeasurements();
}
async function loadStructure(id, keepPose = false) {
  const token = ++state.loadToken;
  if (!id) return;
  setError(""); $("loaded-status").textContent = "Reading coordinates…";
  try {
    const data = await fetchJSON(`/api/structure?id=${encodeURIComponent(id)}`);
    if (token !== state.loadToken) return;
    if (!Array.isArray(data.symbols) || !Array.isArray(data.positions) || !data.symbols.length || data.symbols.length !== data.positions.length || data.positions.some(p => !Array.isArray(p) || p.length !== 3 || !p.every(finite))) throw new Error("The record does not contain a valid, nonempty atomic coordinate set.");
    state.structure = data; state.selected = [];
    if (!keepPose) { state.pose = null; $("endpoint-controls").hidden = true; }
    updateStructureDetails(); viewer.reset(); $("canvas-empty").hidden = true;
    $("loaded-status").textContent = "Coordinates loaded from local evidence";
  } catch (error) {
    if (token !== state.loadToken) return;
    setError(error.message); $("loaded-status").textContent = "Structure could not be loaded";
    if (state.structure) $("structure-select").value = state.structure.id;
    if (!state.structure) $("canvas-empty").hidden = false;
  }
}

class MolecularViewer {
  constructor(canvas) {
    this.canvas = canvas; this.ctx = canvas.getContext("2d"); this.projected = []; this.rx = 1.37; this.ry = -0.45; this.rz = -0.20; this.zoom = 1; this.drag = null;
    new ResizeObserver(() => this.draw()).observe(canvas);
    canvas.addEventListener("pointerdown", e => { if (e.button !== 0) return; this.drag = { x: e.clientX, y: e.clientY, startX: e.clientX, startY: e.clientY, moved: false }; canvas.setPointerCapture(e.pointerId); });
    canvas.addEventListener("pointermove", e => {
      if (this.drag) {
        const dx = e.clientX - this.drag.x, dy = e.clientY - this.drag.y;
        if (Math.hypot(e.clientX - this.drag.startX, e.clientY - this.drag.startY) > 3) this.drag.moved = true;
        this.ry += dx * .009; this.rx += dy * .009; this.drag.x = e.clientX; this.drag.y = e.clientY;
        $("atom-tooltip").hidden = true; this.draw();
      } else {
        const rect = canvas.getBoundingClientRect(), hit = this.hit(e.clientX - rect.left, e.clientY - rect.top);
        const tip = $("atom-tooltip"); tip.hidden = !hit;
        if (hit) { tip.textContent = `${state.structure.symbols[hit.index]} ${hit.index}${atomRole(hit.index) ? " · " + atomRole(hit.index) : ""}`; tip.style.left = `${Math.min(rect.width - 180, Math.max(8, hit.x + 13))}px`; tip.style.top = `${Math.max(44, hit.y - 33)}px`; }
      }
    });
    canvas.addEventListener("pointerup", e => { if (this.drag && !this.drag.moved) { const rect = canvas.getBoundingClientRect(); const hit = this.hit(e.clientX - rect.left, e.clientY - rect.top); if (hit) selectAtom(hit.index); } this.drag = null; });
    canvas.addEventListener("pointercancel", () => { this.drag = null; });
    canvas.addEventListener("pointerleave", () => { $("atom-tooltip").hidden = true; });
    canvas.addEventListener("wheel", e => { e.preventDefault(); this.changeZoom(Math.exp(-e.deltaY * .001)); }, { passive: false });
    canvas.addEventListener("keydown", e => {
      if (e.key === "ArrowLeft") this.ry -= .12; else if (e.key === "ArrowRight") this.ry += .12; else if (e.key === "ArrowUp") this.rx -= .12; else if (e.key === "ArrowDown") this.rx += .12;
      else if (e.key === "+" || e.key === "=") this.changeZoom(1.15); else if (e.key === "-") this.changeZoom(1 / 1.15); else if (e.key === "Home") this.reset(); else return;
      e.preventDefault(); this.draw();
    });
  }
  changeZoom(factor) { this.zoom = Math.max(.35, Math.min(4, this.zoom * factor)); this.draw(); }
  reset() { this.rx = 1.37; this.ry = -.45; this.rz = -.20; this.zoom = 1; this.draw(); }
  rotate(p) {
    const [x, y, z] = p, cx = Math.cos(this.rx), sx = Math.sin(this.rx), cy = Math.cos(this.ry), sy = Math.sin(this.ry), cz = Math.cos(this.rz), sz = Math.sin(this.rz);
    const y1 = y * cx - z * sx, z1 = y * sx + z * cx, x2 = x * cy + z1 * sy, z2 = -x * sy + z1 * cy;
    return [x2 * cz - y1 * sz, x2 * sz + y1 * cz, z2];
  }
  hit(x, y) { return [...this.projected].reverse().find(p => Math.hypot(p.x - x, p.y - y) < p.radius + 5); }
  draw() {
    const canvas = this.canvas, ctx = this.ctx;
    if (!ctx) return;
    const width = canvas.clientWidth, height = canvas.clientHeight;
    if (!width || !height) return;
    const dpr = window.devicePixelRatio || 1;
    if (canvas.width !== Math.round(width * dpr) || canvas.height !== Math.round(height * dpr)) { canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr); }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, width, height);
    const structure = state.structure; this.projected = [];
    if (!structure) return;
    const meta = structureMeta(structure), fixed = new Set(array(structure.fixed_indices)), reactive = new Set(Object.values(transfer()).filter(Number.isInteger));
    const tool = new Set(array(meta.tool_indices)), substrate = new Set(array(meta.substrate_indices));
    const positions = structure.positions, center = [0, 1, 2].map(i => positions.reduce((sum, p) => sum + p[i], 0) / positions.length);
    const relative = positions.map(p => p.map((v, i) => v - center[i]));
    const extent = Math.max(1, ...relative.map(p => Math.hypot(...p))) + .6;
    const scale = Math.min(width, height) * .43 / extent * this.zoom;
    const showH = $("show-hydrogen").checked, group = $("group-select").value;
    const visible = i => (showH || structure.symbols[i] !== "H") && (group === "all" || (group === "tool" ? tool : substrate).has(i));
    const projected = relative.map((p, index) => {
      const [x, y, z] = this.rotate(p); const symbol = structure.symbols[index];
      return { x: width / 2 + x * scale, y: height / 2 + y * scale - 3, z, index, radius: Math.max(3, scale * (symbol === "H" ? .22 : .37)), symbol };
    });
    // A sparse coordinate-plane grid provides orientation, never geometric data.
    ctx.strokeStyle = "#cad4e31f"; ctx.lineWidth = 1;
    for (let i = -4; i <= 4; i++) { const x = width / 2 + i * 40; ctx.beginPath(); ctx.moveTo(x, 55); ctx.lineTo(x, height - 45); ctx.stroke(); }
    for (let i = -3; i <= 3; i++) { const y = height / 2 + i * 40; ctx.beginPath(); ctx.moveTo(40, y); ctx.lineTo(width - 40, y); ctx.stroke(); }
    const items = [];
    array(structure.bonds).forEach(bond => { const a = projected[bond[0]], b = projected[bond[1]]; if (a && b && visible(a.index) && visible(b.index)) items.push({ type: "bond", z: (a.z + b.z) / 2, a, b }); });
    projected.filter(p => visible(p.index)).forEach(p => items.push({ type: "atom", z: p.z, p }));
    items.sort((a, b) => a.z - b.z);
    for (const item of items) {
      if (item.type === "bond") {
        const { a, b } = item;
        ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.lineWidth = Math.max(2, scale * .095); ctx.lineCap = "round";
        ctx.strokeStyle = "#a4afc1"; ctx.stroke();
        ctx.beginPath(); ctx.moveTo(a.x - .7, a.y - .7); ctx.lineTo(b.x - .7, b.y - .7); ctx.lineWidth = Math.max(.7, scale * .025); ctx.strokeStyle = "#e4e9ef"; ctx.stroke();
      } else {
        const p = item.p; const isH = p.symbol === "H", isReactive = reactive.has(p.index), isTool = tool.has(p.index), isSubstrate = substrate.has(p.index);
        if (fixed.has(p.index)) { ctx.beginPath(); ctx.arc(p.x, p.y, p.radius + 4, 0, Math.PI * 2); ctx.strokeStyle = "#17aab4"; ctx.lineWidth = 1.8; ctx.stroke(); }
        if (state.selected.includes(p.index)) { ctx.beginPath(); ctx.arc(p.x, p.y, p.radius + 7, 0, Math.PI * 2); ctx.strokeStyle = "#7b59c2"; ctx.lineWidth = 2; ctx.stroke(); }
        const gradient = ctx.createRadialGradient(p.x - p.radius * .35, p.y - p.radius * .38, p.radius * .04, p.x, p.y, p.radius);
        gradient.addColorStop(0, isReactive ? "#ffedb7" : isH ? "#ffffff" : "#9aa5b9"); gradient.addColorStop(.45, isReactive ? "#e4bd74" : isH ? "#edf0f5" : "#57647b"); gradient.addColorStop(1, isReactive ? "#bd8d41" : isH ? "#c8d0de" : "#2e394e");
        ctx.fillStyle = gradient; ctx.beginPath(); ctx.arc(p.x, p.y, p.radius, 0, Math.PI * 2); ctx.fill();
        ctx.strokeStyle = isReactive ? "#bc903f" : isTool ? "#8c80b5" : isSubstrate ? "#6f9ea7" : "#8f9aac"; ctx.lineWidth = isH ? .6 : 1.2; ctx.stroke();
      }
    }
    this.projected = projected.filter(p => visible(p.index)).sort((a, b) => a.z - b.z);
    const showLabels = $("show-indices").checked;
    for (const p of this.projected) {
      if (!showLabels && !state.selected.includes(p.index)) continue;
      ctx.font = `600 ${showLabels ? 9 : 11}px -apple-system, sans-serif`; ctx.textAlign = "center"; ctx.textBaseline = "middle";
      if (showLabels) { ctx.fillStyle = p.symbol === "H" || reactive.has(p.index) ? "#594d38" : "#ffffff"; ctx.fillText(String(p.index), p.x, p.y); }
      else { const text = `${p.symbol} ${p.index}`; const tw = ctx.measureText(text).width; ctx.fillStyle = "#ffffffed"; ctx.fillRect(p.x - tw / 2 - 5, p.y - p.radius - 24, tw + 10, 17); ctx.fillStyle = "#7153af"; ctx.fillText(text, p.x, p.y - p.radius - 15); }
    }
    if (state.selected.length === 2) {
      const [a, b] = state.selected.map(i => projected[i]);
      if (a && b && visible(a.index) && visible(b.index)) {
        ctx.beginPath(); ctx.setLineDash([4, 4]); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.strokeStyle = "#8966be"; ctx.lineWidth = 1; ctx.stroke(); ctx.setLineDash([]);
      }
    }
    const groupLabels = [[tool, "TOOL", "#8983ae"], [substrate, "SUBSTRATE", "#759aa5"]].map(([indices, label, color]) => {
      const points = projected.filter(p => indices.has(p.index) && visible(p.index));
      if (!points.length) return null;
      return { label, color, count: points.length, y: Math.max(80, Math.min(height - 72, points.reduce((sum, p) => sum + p.y, 0) / points.length)) };
    }).filter(Boolean);
    if (groupLabels.length === 2 && Math.abs(groupLabels[0].y - groupLabels[1].y) < 35) { groupLabels[0].y -= 20; groupLabels[1].y += 20; }
    groupLabels.forEach(item => { ctx.textAlign = "left"; ctx.textBaseline = "middle"; ctx.fillStyle = item.color; ctx.font = "600 8px -apple-system, sans-serif"; ctx.fillText(item.label, 18, item.y); ctx.font = "8px -apple-system, sans-serif"; ctx.fillText(`${item.count} visible atoms`, 18, item.y + 13); });
    this.drawAxes(ctx, width, height);
  }
  drawAxes(ctx, width, height) {
    const origin = [width - 43, height - 60];
    [["x", "#b89191", [1, 0, 0]], ["y", "#8da797", [0, 1, 0]], ["z", "#899bb9", [0, 0, 1]]].forEach(([label, color, axis]) => {
      const v = this.rotate(axis); ctx.beginPath(); ctx.moveTo(...origin); ctx.lineTo(origin[0] + v[0] * 20, origin[1] + v[1] * 20); ctx.strokeStyle = color; ctx.lineWidth = 1.5; ctx.stroke(); ctx.font = "8px sans-serif"; ctx.fillStyle = color; ctx.textAlign = "center"; ctx.fillText(label, origin[0] + v[0] * 26, origin[1] + v[1] * 26);
    });
  }
}
const viewer = new MolecularViewer($("molecule-canvas"));

function switchTab(tab) {
  state.tab = tab;
  const labels = { structure: ["Atomic structure", "Inspect the candidate, one atom at a time."], poses: ["Pose study", "Compare recorded starting geometries and their calculation status."], calculations: ["Calculation evidence", "Executed results, numerical diagnostics, and remaining validation."], references: ["Reference comparison", "Examine method differences and the sensitivity to electronic solutions."], path: ["Path & vibrations", "Reaction connectivity and stationary-point evidence."] };
  $("page-title").textContent = labels[tab][0]; $("page-description").textContent = labels[tab][1];
  document.querySelectorAll(".nav-button").forEach(button => { const active = button.dataset.tab === tab; button.classList.toggle("active", active); if (active) button.setAttribute("aria-current", "page"); else button.removeAttribute("aria-current"); });
  $("molecule-section").hidden = !["structure", "poses"].includes(tab);
  $("structure-summary").hidden = !["structure", "poses"].includes(tab);
  ["poses", "calculations", "references", "path"].forEach(name => { $(`${name}-section`).hidden = name !== tab; });
  if (tab === "poses" && !state.pose && state.campaign?.entries?.length) selectPose(state.campaign.entries.find(entry => entry.initial_structure_id === state.structure?.id) || state.campaign.entries[4] || state.campaign.entries[0]);
  viewer.draw();
}
function selectPose(entry, endpoint = "initial") {
  state.pose = entry; state.endpoint = endpoint;
  $("endpoint-controls").hidden = false;
  $("initial-button").classList.toggle("selected", endpoint === "initial"); $("final-button").classList.toggle("selected", endpoint === "final");
  $("initial-button").disabled = !entry.initial_structure_id; $("final-button").disabled = !entry.final_structure_id;
  const id = entry[`${endpoint}_structure_id`];
  if (id) loadStructure(id, true); else setError(`No ${endpoint} structure is available for this pose.`);
  document.querySelectorAll(".pose-row").forEach(row => row.classList.toggle("selected-row", row.dataset.pose === entry.id));
}
function renderPoses() {
  const section = $("poses-section"); section.replaceChildren();
  const campaigns = array(state.catalog.campaigns);
  $("pose-count").textContent = campaigns.reduce((sum, item) => sum + array(item.entries).length, 0);
  if (!campaigns.length) { section.append(empty("No campaign loaded", "Start the workbench with a configured local campaign directory to inspect a pose study.")); return; }
  if (!state.campaign) state.campaign = campaigns[0];
  const card = node("div", "content-card"), heading = node("div", "content-heading"), intro = node("div");
  append(intro, node("h2", "", "Pose campaign"), node("p", "", "Select a row to inspect its saved endpoints. Pending poses have no calculated result."));
  const select = node("select"); select.setAttribute("aria-label", "Campaign");
  campaigns.forEach(c => select.add(new Option(c.label || c.id, c.id))); select.value = state.campaign.id;
  select.addEventListener("change", () => { state.campaign = campaigns.find(c => c.id === select.value); state.pose = null; renderPoses(); if (state.campaign.entries?.length) selectPose(state.campaign.entries[0]); });
  append(heading, intro, select); card.append(heading);
  if (state.campaign.load_error) card.append(node("div", "error-banner", state.campaign.load_error));
  array(state.campaign.notices).forEach(value => card.append(node("div", "notice-box", value)));
  array(state.campaign.integrity_errors).forEach(value => card.append(node("div", "error-banner", value)));
  const table = node("table"), head = node("thead"), headerRow = node("tr");
  ["POSE", "SEPARATION · Å", "OFFSET · Å", "STATUS", "ENDPOINTS"].forEach(t => headerRow.append(node("th", "", t))); head.append(headerRow); table.append(head);
  const body = node("tbody");
  array(state.campaign.entries).forEach((entry, index) => {
    const row = node("tr", "pose-row"); row.dataset.pose = entry.id; row.tabIndex = 0; row.setAttribute("aria-label", `Inspect ${entry.label || entry.id}`);
    const pose = object(entry.pose); const sep = pose.separation_angstrom ?? entry.separation_angstrom; const offset = pose.offset_angstrom ?? pose.lateral_offset_angstrom ?? entry.offset_angstrom;
    append(row, node("td", "", entry.label || entry.id || `Pose ${index + 1}`), node("td", "", fmt(sep, 2)), node("td", "", signed(offset, 2)), append(node("td"), statusPill(entry.status)), node("td", "", `${entry.initial_structure_id ? "Initial" : "—"}${entry.final_structure_id ? " / final guess" : ""}`));
    if (entry.load_error || entry.evidence_error) row.cells[3].append(node("div", "small", entry.load_error || entry.evidence_error));
    array(entry.integrity_errors).forEach(value => row.cells[3].append(node("div", "small", value)));
    if (entry.status_verified === false) row.cells[3].append(node("div", "small", "Provenance unverified"));
    row.addEventListener("click", () => selectPose(entry)); row.addEventListener("keydown", e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); selectPose(entry); } }); body.append(row);
  });
  table.append(body); card.append(append(node("div", "table-wrap"), table));
  if (!array(state.campaign.entries).length) card.append(node("p", "muted small", "No readable pose entries in this campaign."));
  append(card, node("p", "muted small", "Separation and lateral offset are input settings. The actual donor–apex distance is measured from the selected coordinates above."), sourceLink(state.campaign));
  section.append(card);
}
function calculationCard(record) {
  const card = node("article", "record-card"), values = object(record.structure), q = object(values.quantum_diagnostics), settings = object(record.quantum_settings || q.settings);
  const header = node("div", "record-top"); append(header, append(node("div"), node("h3", "", record.label || record.id), node("p", "", `${human(record.stage)} · ${human(record.state)} geometry`)), statusPill(record.status)); card.append(header);
  card.append(append(node("div", "energy-value", fmt(values.energy_ev, 6)), node("small", "", " eV")), node("div", "metric-label", "Recorded total electronic energy"));
  const metrics = node("div", "metric-grid");
  append(metrics, metric("Maximum free force · eV/Å", fmt(values.free_force_max_ev_per_angstrom, 6)), metric("Elapsed time · seconds", fmt(record.elapsed_seconds ?? q.elapsed_seconds, 1)), metric("Determinant ⟨S²⟩", fmt(q.s2, 6)), metric("Expected S(S+1)", fmt(q.expected_s2, 3))); card.append(metrics);
  const method = `${typeof settings.xc === "string" ? settings.xc.toUpperCase() : "Not recorded"} / ${settings.basis || "—"}${settings.dispersion ? " · " + settings.dispersion : ""}`;
  card.append(detailList([["Method", method], ["DFT SCF initial guess", settings.scf_initial_guess ?? "Not recorded"], ["Density fitting", yesNo(settings.density_fit)], ["SCF converged", yesNo(q.scf_converged)], ["Gradient completed", yesNo(q.gradient_completed)], ["Stability checked", yesNo(q.stability_checked)], ["SCF cycles", q.scf_cycles], ["Basis functions", q.basis_functions], ["Threads requested / effective", `${q.requested_pyscf_threads ?? settings.threads ?? "—"} / ${q.effective_pyscf_threads ?? "—"}`], ["Grid level", settings.grid_level], ["SCF tolerance", finite(settings.conv_tol) ? settings.conv_tol.toExponential() : "Not recorded"], ["Charge / spin (2S)", `${settings.charge ?? "—"} / ${settings.spin ?? "—"}`]]));
  const groups = object(values.anchor_force_groups);
  if (array(values.anchor_force_sum_ev_per_angstrom).length === 3) {
    const sums = detailList([["Global internal anchor sum · eV/Å", `[${values.anchor_force_sum_ev_per_angstrom.map(v => fmt(v, 4)).join(", ")}]`]]);
    sums.querySelectorAll("dd").forEach(el => el.classList.add("anchor-vector")); card.append(sums);
  }
  if (Object.keys(groups).length) {
    const list = detailList(Object.entries(groups).map(([name, value]) => [`${human(name)} holding · eV/Å`, array(value.external_holding_force_ev_per_angstrom).length === 3 ? `[${value.external_holding_force_ev_per_angstrom.map(v => fmt(v, 4)).join(", ")}]` : "Not recorded"]));
    list.querySelectorAll("dd").forEach(el => el.classList.add("anchor-vector")); card.append(list);
    card.append(node("p", "muted small", "Holding-force vectors are the negatives of the group internal-force sums. Opposing loads can cancel in a global sum."));
  } else card.append(node("p", "muted small", "Separate group anchor loads were not recorded. The global sum can hide opposing loads."));
  const interpretation = object(record.interpretation);
  if (interpretation.original_run_code_context || interpretation.timing_caveat || interpretation.event_log_status) {
    const details = node("details"); details.append(node("summary", "", "Archive and timing context"));
    [interpretation.original_run_code_context, interpretation.event_log_status, interpretation.timing_caveat].filter(Boolean).forEach(text => details.append(node("p", "muted small", text)));
    append(details, sourceLink({ source_href: record.interpretation_href }, "View archive interpretation ↗")); card.append(details);
  }
  const validation = object(record.validation);
  card.append(node("p", "notice-box", validation.accuracy_claim || "This calculation does not establish a relaxed design, reaction barrier, or validated tool."));
  const missing = array(validation.missing_evidence);
  if (missing.length) { const details = node("details"), list = node("ul", "list-note"); details.append(node("summary", "", `Remaining evidence (${missing.length})`)); missing.forEach(value => list.append(node("li", "", value))); append(details, list); card.append(details); }
  if (record.error) card.append(node("p", "error-banner", typeof record.error === "string" ? record.error : JSON.stringify(record.error)));
  if (record.load_error) card.append(node("p", "error-banner", record.load_error));
  const diagnostic = node("details"); append(diagnostic, node("summary", "", "Recorded settings and diagnostics"), node("pre", "", JSON.stringify({ quantum_settings: settings, quantum_diagnostics: q, validation }, null, 2))); card.append(diagnostic);
  append(card, sourceLink(record));
  if (record.structure_id) { const button = node("button", "text-button", "Inspect this structure →"); button.style.marginLeft = "18px"; button.addEventListener("click", () => { switchTab("structure"); loadStructure(record.structure_id); }); card.append(button); }
  return card;
}
function renderCalculations() {
  const section = $("calculations-section"); section.replaceChildren();
  const records = array(state.catalog.calculations);
  if (!records.length) { section.append(empty("No calculation records", "No result records are available in the configured evidence root.")); return; }
  section.append(node("div", "notice-box", "Completed is a calculation status. Single-point energies are not reaction barriers or relaxed designs. Direct and density-fitting runs may overlap in time; these elapsed times are not a controlled speed benchmark."));
  const grid = node("div", "record-grid"); records.forEach(record => grid.append(calculationCard(record))); section.append(grid);
}
function referenceCard(record) {
  const card = node("article", "record-card reference-card"), cc = object(record.cc_settings), computed = object(record.computed), dft = object(computed.dft), coupled = object(computed.ccsd_t);
  const annotation = object(record.retrospective_annotation), annotationBound = record.retrospective_annotation_status === "verified source binding";
  const inferredCCGuess = annotationBound && typeof annotation.inferred_cc_initial_guess === "string" ? annotation.inferred_cc_initial_guess : null;
  const ccGuessLabel = cc.scf_initial_guess ?? (inferredCCGuess ? `${inferredCCGuess} (inferred; not originally recorded)` : "Not recorded");
  const heading = node("div", "record-top"); append(heading, append(node("div"), node("h3", "", record.label || record.id), node("p", "", `CC reference initial guess: ${ccGuessLabel} · ${record.basis || cc.basis || "basis not recorded"}`)), statusPill(record.status)); card.append(heading);
  append(card, append(node("div", "energy-value", signed(coupled.nominal_ts_relative_energy_kcal_per_mol, 4)), node("small", "", " kcal/mol")), node("p", "computed-label", "CCSD(T): source transition geometry relative to separated reactants"));
  card.append(append(node("div", "metric-grid"), metric("DFT relative energy · kcal/mol", signed(dft.nominal_ts_relative_energy_kcal_per_mol, 4)), metric("DFT − CCSD(T) · kcal/mol", signed(computed.dft_minus_ccsd_t?.nominal_ts_relative_energy_kcal_per_mol, 4))));
  const table = node("table"), head = node("tr"); ["SPECIES", "HF ⟨S²⟩", "KS ⟨S²⟩"].forEach(label => head.append(node("th", "", label))); table.append(append(node("thead"), head));
  const body = node("tbody");
  Object.entries(object(record.species)).forEach(([name, item]) => { append(body, append(node("tr"), node("td", "", human(name)), node("td", "", fmt(item.cc_hf_s2, 4)), node("td", "", fmt(item.dft_s2, 4)))); }); table.append(body); card.append(append(node("div", "table-wrap"), table));
  card.append(node("p", "muted small", "Spin values describe the HF and Kohn–Sham reference determinants. They do not establish correlated CC spin purity or matching electronic states."));
  card.append(detailList([["Recorded CC SCF initial guess", cc.scf_initial_guess ?? "Not recorded"], ["Recorded DFT SCF initial guess", object(record.dft_settings).scf_initial_guess ?? "Not recorded"], ["Geometry/basis pairing verified", yesNo(record.geometry_basis_pairing_verified)], ["Electronic-state identity verified", yesNo(record.electronic_state_identity_verified)], ["Ground state verified", yesNo(record.ground_state_verified)], ["First-order saddle verified", yesNo(record.saddle_verified)], ["Method accuracy validated", yesNo(record.method_accuracy_validated)]]));
  if (record.interpretation) append(card, node("div", "metric-label", "Original record interpretation"), node("p", "muted small", record.interpretation));
  const limitations = array(record.limitations);
  if (limitations.length) { const details = node("details"), list = node("ul", "list-note"); limitations.forEach(value => list.append(node("li", "", value))); append(details, node("summary", "", "Interpretation and limitations"), list); card.append(details); }
  if (Object.keys(annotation).length || record.retrospective_annotation_status === "invalid source binding" || array(record.annotation_errors).length) {
    const section = node("section", "historical-annotation"); section.setAttribute("aria-label", "Retrospective reference annotation");
    append(section, node("h3", "", "Retrospective annotation"), node("div", "metric-label", `${annotation.annotation_date || "Date not recorded"} · ${record.retrospective_annotation_status || "Source binding not recorded"}`));
    if (annotationBound) {
      if (!cc.scf_initial_guess && inferredCCGuess) section.append(node("p", "small", `Inferred CC initial guess: ${inferredCCGuess}. The original record did not serialize this choice.`));
      if (annotation.inference_basis) section.append(node("p", "muted small", annotation.inference_basis));
      if (annotation.interpretation) section.append(node("p", "small", annotation.interpretation));
      section.append(detailList([["Electronic-state identity verified", yesNo(annotation.electronic_state_identity_verified)], ["Ground state verified", yesNo(annotation.ground_state_verified)], ["First-order saddle verified", yesNo(annotation.saddle_verified)], ["Method accuracy validated", yesNo(annotation.method_accuracy_validated)], ["Original record modified", yesNo(annotation.original_record_modified)]]));
      if (annotation.historical_data_policy) section.append(node("p", "muted small", annotation.historical_data_policy));
    } else section.append(node("p", "error-banner", "This annotation is not bound to the loaded original record. Its scientific interpretation and inferred settings are not accepted as metadata."));
    array(record.annotation_errors).forEach(message => section.append(node("p", "error-banner", message)));
    append(section, sourceLink({ source_href: record.retrospective_annotation_href }, "View retrospective annotation source ↗")); card.append(section);
  }
  if (record.load_error) card.append(node("p", "error-banner", record.load_error));
  append(card, sourceLink(record)); return card;
}
function renderReferences() {
  const section = $("references-section"); section.replaceChildren(); const records = array(state.catalog.references);
  if (!records.length) { section.append(empty("No reference comparisons", "Paired method records have not been loaded.")); return; }
  section.append(node("div", "notice-box", "These are fixed-geometry relative energies, not verified reaction barriers. Different initial guesses can converge to different electronic solutions. Numerical convergence does not identify the ground state."));
  const grid = node("div", "record-grid"); records.forEach(record => grid.append(referenceCard(record))); section.append(grid);
}
function renderPath() {
  const section = $("path-section"); section.replaceChildren();
  section.append(empty("Not computed for this candidate", "No connected reaction path or stationary-point vibration record is loaded for the 53-atom candidate. Endpoint guesses and single-point forces cannot establish a transition state or barrier.", "⌁"));
  const card = node("div", "content-card"); append(card, node("h2", "", "What would establish this evidence?"), node("p", "reference-lead", "A recorded path, converged force residuals, a verified first-order saddle and displacement-step checks, followed by connectivity to the intended reactant and product. Electronic-state and method validation remain separate requirements.")); section.append(card);
}
async function initialize() {
  try {
    state.catalog = await fetchJSON("/api/catalog");
    const select = $("structure-select"); select.replaceChildren();
    array(state.catalog.structures).forEach(s => select.add(new Option(s.label || s.id, s.id)));
    const notices = array(state.catalog.notices).map(value => typeof value === "string" ? value : value.message || JSON.stringify(value));
    if (notices.length) { $("notice-banner").textContent = notices.join(" · "); $("notice-banner").hidden = false; }
    renderPoses(); renderCalculations(); renderReferences(); renderPath();
    const id = state.catalog.default_structure_id || state.catalog.structures?.[0]?.id;
    if (id) await loadStructure(id); else { $("canvas-empty").hidden = false; $("loaded-status").textContent = "No structures loaded"; $("structure-label").textContent = "No structures available"; $("formula").textContent = "—"; }
  } catch (error) { setError(`Cannot load the local workbench: ${error.message}`); $("loaded-status").textContent = "Local evidence unavailable"; $("canvas-empty").hidden = false; }
}
document.querySelectorAll(".nav-button").forEach(button => button.addEventListener("click", () => switchTab(button.dataset.tab)));
$("structure-select").addEventListener("change", e => loadStructure(e.target.value));
$("reset-view").addEventListener("click", () => viewer.reset());
$("zoom-in").addEventListener("click", () => viewer.changeZoom(1.2)); $("zoom-out").addEventListener("click", () => viewer.changeZoom(1 / 1.2));
["show-indices", "show-hydrogen", "group-select"].forEach(id => $(id).addEventListener("change", () => viewer.draw()));
$("clear-selection").addEventListener("click", () => { state.selected = []; updateMeasurements(); });
["atom-a", "atom-b"].forEach(id => $(id).addEventListener("change", () => { const a = $("atom-a").value, b = $("atom-b").value; state.selected = [a, b].filter(value => value !== "").map(Number); if (state.selected[0] === state.selected[1]) state.selected.pop(); updateMeasurements(); }));
$("initial-button").addEventListener("click", () => { if (state.pose) selectPose(state.pose, "initial"); });
$("final-button").addEventListener("click", () => { if (state.pose) selectPose(state.pose, "final"); });
initialize();
