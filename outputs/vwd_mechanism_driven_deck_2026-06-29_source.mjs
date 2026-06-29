import fs from "node:fs/promises";
import path from "node:path";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const OUT_DIR = "/Users/lucachangretta/Agent_Workspace/VWF-ETHos/outputs";
const PREVIEW_DIR = "/tmp/codex-presentations/vwd_mechanism_deck_20260629/tmp/preview";
const LAYOUT_DIR = "/tmp/codex-presentations/vwd_mechanism_deck_20260629/tmp/layout";
const FINAL = path.join(OUT_DIR, "vwd_mechanism_driven_agentic_diagnosis_2026-06-29.pptx");

const W = 1280;
const H = 720;
const C = {
  ink: "#111111",
  muted: "#555555",
  light: "#F1F2F4",
  mid: "#D9DDE3",
  line: "#B8BCC4",
  orange: "#FF6B35",
  blue: "#2563EB",
  green: "#168A5A",
  red: "#C2410C",
  purple: "#6D28D9",
};

function pos(left, top, width, height) {
  return { left, top, width, height };
}

function addText(slide, text, p, opts = {}) {
  const shape = slide.shapes.add({
    geometry: "textbox",
    position: p,
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  shape.text = text;
  shape.text.style = {
    fontSize: opts.size ?? 22,
    bold: opts.bold ?? false,
    color: opts.color ?? C.ink,
    alignment: opts.align ?? "left",
  };
  return shape;
}

function addBox(slide, p, opts = {}) {
  return slide.shapes.add({
    geometry: opts.geometry ?? "rect",
    position: p,
    fill: opts.fill ?? C.light,
    line: { style: "solid", fill: opts.line ?? C.mid, width: opts.lineWidth ?? 1 },
  });
}

function title(slide, t, sub = "") {
  addText(slide, t, pos(48, 38, 940, 74), { size: 38, bold: true });
  if (sub) addText(slide, sub, pos(50, 112, 920, 30), { size: 17, color: C.muted });
  addText(slide, "VWF-ETHos · 2026-06-29", pos(1035, 48, 190, 24), { size: 14, color: C.muted, align: "right" });
}

function foot(slide, txt) {
  addText(slide, txt, pos(48, 684, 1120, 20), { size: 11, color: "#777777" });
}

function metric(slide, label, value, p, color = C.ink, sub = "") {
  addText(slide, value, pos(p.left, p.top, p.width, 58), { size: 42, bold: true, color });
  addText(slide, label, pos(p.left, p.top + 58, p.width, 28), { size: 16, color: C.muted });
  if (sub) addText(slide, sub, pos(p.left, p.top + 86, p.width, 40), { size: 14, color: C.muted });
}

function bulletList(slide, items, p, opts = {}) {
  const text = items.map((x) => `• ${x}`).join("\n");
  return addText(slide, text, p, { size: opts.size ?? 20, color: opts.color ?? C.ink });
}

function smallTable(slide, headers, rows, p, widths, opts = {}) {
  const rowH = opts.rowH ?? 38;
  const headH = opts.headH ?? 34;
  let x = p.left;
  addBox(slide, pos(p.left, p.top, p.width, headH), { fill: C.ink, line: C.ink });
  headers.forEach((h, i) => {
    addText(slide, h, pos(x + 8, p.top + 8, widths[i] - 16, headH - 8), { size: opts.headerSize ?? 14, bold: true, color: "#FFFFFF" });
    x += widths[i];
  });
  rows.forEach((r, ri) => {
    const y = p.top + headH + ri * rowH;
    addBox(slide, pos(p.left, y, p.width, rowH), { fill: ri % 2 ? "#FFFFFF" : C.light, line: C.mid });
    let cx = p.left;
    r.forEach((cell, ci) => {
      const color = opts.cellColor?.(ri, ci, cell) ?? C.ink;
      addText(slide, String(cell), pos(cx + 8, y + 9, widths[ci] - 16, rowH - 8), { size: opts.size ?? 14, color, bold: opts.boldCell?.(ri, ci, cell) ?? false });
      cx += widths[ci];
    });
  });
}

function arrow(slide, from, to, label = "") {
  slide.shapes.add({
    geometry: "line",
    position: { left: from.x, top: from.y, width: to.x - from.x, height: to.y - from.y },
    line: { style: "solid", fill: C.line, width: 2, beginArrowType: "none", endArrowType: "triangle" },
  });
  if (label) addText(slide, label, pos((from.x + to.x) / 2 - 50, (from.y + to.y) / 2 - 22, 100, 22), { size: 12, color: C.muted, align: "center" });
}

async function writeBlob(file, blob) {
  await fs.writeFile(file, new Uint8Array(await blob.arrayBuffer()));
}

async function main() {
  await fs.mkdir(OUT_DIR, { recursive: true });
  await fs.mkdir(PREVIEW_DIR, { recursive: true });
  await fs.mkdir(LAYOUT_DIR, { recursive: true });

  const deck = Presentation.create({ slideSize: { width: W, height: H } });

  // 1
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    addText(s, "以机制驱动的 VWD", pos(48, 72, 860, 86), { size: 58, bold: true });
    addText(s, "智能体辅助分型/诊断系统", pos(48, 154, 1010, 86), { size: 58, bold: true });
    addText(s, "研究原型 · 面向 VWF 变异的可解释 subtype 证据融合", pos(52, 266, 900, 38), { size: 23, color: C.muted });
    addBox(s, pos(48, 390, 1184, 2), { fill: C.ink, line: C.ink });
    metric(s, "最新实际 Eval", "340", pos(52, 440, 210, 120), C.blue, "Type 1 + Type 2");
    metric(s, "Type 2 样本", "225", pos(310, 440, 210, 120), C.ink, "2A/2B/2M/2N");
    metric(s, "Type 2 recall", "54.2%", pos(568, 440, 240, 120), C.orange, "no hotspot prior");
    metric(s, "下一阶段目标", "机制补全", pos(864, 440, 280, 120), C.green, "Type 1/3 + 2M/2A");
    foot(s, "Method name suggestion: mechanism-driven VWD agentic decision-support system; diagnosis remains clinician-confirmed.");
  }

  // 2
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "为什么叫“机制驱动”", "VWD 的临床 subtype 本身就是病理机制标签；系统目标不是黑箱分类，而是把证据按机制归位。");
    const rows = [
      ["Type 1", "部分定量缺乏", "VWF 水平/分泌/清除"],
      ["Type 3", "几乎完全缺乏", "null / biallelic / severe deficiency"],
      ["2A", "高分子多聚体减少", "A2 稳定性、ADAMTS13、多聚化/分泌"],
      ["2B", "A1-GPIb gain-of-function", "AIM 自抑制松开且结合面保留"],
      ["2M", "功能下降但多聚体保留", "A1-GPIb LOF 或胶原结合 LOF"],
      ["2N", "FVIII 结合下降", "D'/D3-FVIII binding"],
    ];
    smallTable(s, ["分型", "临床/病理定义", "对应机制轴"], rows, pos(58, 172, 720, 378), [120, 240, 360], { rowH: 54, size: 17 });
    addBox(s, pos(840, 188, 330, 310), { fill: "#FFF6F1", line: "#FFD3C2" });
    addText(s, "核心优点", pos(870, 218, 260, 34), { size: 25, bold: true, color: C.orange });
    bulletList(s, [
      "知道缺什么：每个 subtype 可追溯到具体功能轴",
      "知道错在哪：错误 case 能落到机制冲突或证据缺失",
      "知道怎么补：Boltz/AF3/MD/SMD 分别服务不同机制",
      "避免只看总 accuracy 的假安全感",
    ], pos(870, 274, 260, 170), { size: 18 });
    foot(s, "Source basis: ASH/ISTH/NHF/WFH diagnosis guideline; Haematologica 2026 classification review.");
  }

  // 3
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "最新 Eval 直接结果", "不混历史数据：这里仅展示当前 pipeline 实际进入评估的 clean supported labels。");
    smallTable(s, ["真实标签", "n", "当前状态"], [
      ["1", "115", "进入 Eval；但 Type 1 规则尚不完整"],
      ["2A", "118", "进入 Eval"],
      ["2B", "38", "进入 Eval"],
      ["2M", "53", "进入 Eval"],
      ["2N", "16", "进入 Eval"],
      ["3", "44", "干净库存存在；当前分类器未纳入"],
    ], pos(58, 168, 510, 358), [130, 90, 290], { rowH: 54, size: 17 });
    deck.slides.items.at(-1).charts.add("bar", {
      position: pos(640, 168, 540, 330),
      categories: ["1", "2A", "2B", "2M", "2N", "3*"],
      series: [{ name: "样本数", values: [115, 118, 38, 53, 16, 44], fill: C.blue }],
      hasLegend: false,
      dataLabels: { showValue: true, position: "outEnd" },
      yAxis: { majorGridlines: { style: "solid", fill: C.mid, width: 1 } },
    });
    addText(s, "* Type 3 尚未进入本轮 recall/precision 计算", pos(685, 520, 440, 24), { size: 15, color: C.muted });
    foot(s, "Repo artifact: output/type2m_lof_md_fast_validation_2026-06-29/eval_with_fast_md/eval_v2_summary.csv");
  }

  // 4
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "Type 2: 去掉 hotspot 后，系统更干净但更保守", "with MD + no hotspot prior 后，2B recall 下降，2M recall 上升；这暴露了真正需要动态 GOF 证据的缺口。");
    s.charts.add("bar", {
      position: pos(70, 160, 520, 330),
      categories: ["2A", "2B", "2M", "2N"],
      series: [
        { name: "Recall", values: [59.3, 39.5, 47.2, 75.0], fill: C.orange },
        { name: "Precision", values: [93.3, 88.2, 89.3, 30.0], fill: C.blue },
      ],
      hasLegend: true,
      dataLabels: { showValue: true, position: "outEnd" },
      yAxis: { min: 0, max: 100, majorGridlines: { style: "solid", fill: C.mid, width: 1 } },
    });
    metric(s, "Type 2 overall recall", "122/225", pos(680, 155, 260, 120), C.orange, "54.2%, uncertain=65");
    metric(s, "2M precision", "89.3%", pos(960, 155, 220, 120), C.green, "判出来仍较干净");
    addBox(s, pos(670, 338, 490, 130), { fill: C.light, line: C.mid });
    bulletList(s, [
      "2B: recall 明显下降 = 之前确实依赖位置先验",
      "2M: recall 小幅上升 = G1324R/A1377V 不再被位置先验截走",
      "2N: recall 还行，但 precision 低 = 很多 2A 被吸到 2N",
      "no-hotspot 版本更适合汇报和投稿",
    ], pos(695, 365, 440, 84), { size: 17 });
    foot(s, "Precision from current Type 2 confusion matrix: predicted-class purity; recall from true-label recovery.");
  }

  // 5
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "错误结构：哪里漏、哪里误伤", "这页用于解释为什么报告 recall/precision，而不是只给一个总体准确率。");
    smallTable(s, ["true \\ pred", "2A", "2B", "2M", "2N", "uncertain"], [
      ["2A", "70", "2", "0", "28", "18"],
      ["2B", "1", "15", "3", "0", "19"],
      ["2M", "0", "0", "25", "0", "28"],
      ["2N", "4", "0", "0", "12", "0"],
    ], pos(56, 170, 620, 238), [140, 90, 90, 90, 90, 120], { rowH: 48, size: 18,
      cellColor: (ri, ci, cell) => (ci === ri + 1 ? C.green : (cell === "24" || cell === "28" ? C.red : C.ink)),
      boldCell: (ri, ci) => ci === ri + 1,
    });
    addText(s, "三个直接结论", pos(740, 168, 330, 34), { size: 28, bold: true });
    bulletList(s, [
      "2M 最大问题仍是漏：28/53 仍 uncertain",
      "2B 最大问题变成保守：19/38 uncertain，3/38 被误判 2M",
      "2A 最大问题是被 2N 规则吸走：28 例 2A→2N",
    ], pos(742, 222, 410, 140), { size: 19 });
    addBox(s, pos(738, 410, 390, 98), { fill: "#FFF6F1", line: "#FFD3C2" });
    addText(s, "汇报口径", pos(760, 432, 320, 28), { size: 20, bold: true, color: C.orange });
    addText(s, "当前模型应定位为 Type 2 机制分型原型；Type 1/3 补齐前不报告全分型 global accuracy。", pos(760, 466, 330, 38), { size: 16, color: C.ink });
    foot(s, "Repo artifact: eval_v2_type2_confusion_with_md.csv");
  }

  // 6
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "Case 正例：L1276P 从 uncertain 被 MD 救回 2M", "这页把每个指标、计算方式和 MD 补足的盲区讲清楚。");
    const xs = [54, 278, 502, 726, 950];
    const labs = ["输入", "Boltz 静态结构", "规则融合", "Equilibrium MD", "输出"];
    labs.forEach((lab, i) => {
      addBox(s, pos(xs[i], 158, 180, 106), { fill: i === 3 ? "#FFF6F1" : C.light, line: i === 3 ? "#FFD3C2" : C.mid });
      addText(s, lab, pos(xs[i] + 12, 172, 150, 22), { size: 17, bold: true, color: i === 3 ? C.orange : C.ink });
    });
    arrow(s, { x: 234, y: 211 }, { x: 276, y: 211 });
    arrow(s, { x: 458, y: 211 }, { x: 500, y: 211 });
    arrow(s, { x: 682, y: 211 }, { x: 724, y: 211 });
    arrow(s, { x: 906, y: 211 }, { x: 948, y: 211 });
    addText(s, "aa_change=L1276P\nA1 domain, pos=1276\nno position prior", pos(66, 198, 154, 58), { size: 13 });
    addText(s, "fb z=+0.25\nheparan z=+1.16\naim static=-0.50", pos(290, 198, 154, 58), { size: 13 });
    addText(s, "no MD:\nA1 方向证据弱\n→ uncertain", pos(514, 198, 154, 58), { size: 13 });
    addText(s, "AIM-A1 contacts\n0-5ns: 146.8\n40-50ns: 112.8\nloss=34.0", pos(738, 190, 154, 72), { size: 12 });
    addText(s, "pred=2M\nconf=0.55\nMD LOF evidence", pos(962, 198, 154, 58), { size: 13, color: C.green });

    smallTable(s, ["指标", "怎么算", "本 case 含义"], [
      ["fb_binding_zscore", "Boltz A1+GPIb complex primary score 的变体-同轴分布 z-score", "+0.25：不支持明显 GPIb 结合丢失"],
      ["heparan_zscore", "Boltz A1+heparan model 同法转 z-score", "+1.16：静态表面轴不支持 LOF"],
      ["position prior", "已从分类器移除；不再用 recurrent residue 直接判 2B", "避免标签记忆/审稿质疑"],
      ["md_face_destab", "(0-5ns contacts − 40-50ns contacts) / 20", "1.701：tail 接触显著下降，支持 2M LOF"],
    ], pos(58, 322, 1110, 208), [190, 480, 440], { rowH: 43, size: 13 });
    addText(s, "MD 补的盲区：Boltz 给的是静态复合体置信/结合 proxy，不能看到 closed-state A1 结合面在水环境中是否随时间塌陷；equilibrium MD 用 tail window 观察接触网络保留，给出动态稳定性证据。", pos(86, 552, 1030, 42), { size: 17, bold: true, color: C.ink });
    foot(s, "Position/recurrent-residue prior removed after leakage-risk audit. Sources for mechanism context only: ASH/ISTH/NHF/WFH 2021; PDB 7A6O/Nat Commun 2021.");
  }

  // 7
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "Case 负例：V1316M 显示 equilibrium MD 的歧义", "去掉 hotspot 后，不再位置记忆；但单一 equilibrium MD LOF 轴会误伤部分真实 2B。");
    smallTable(s, ["字段", "值", "解释"], [
      ["输入", "V1316M, A1", "true label = 2B"],
      ["位置先验", "已移除", "不再因为 recurrent residue 直接判 2B"],
      ["Boltz 轴", "fb=-0.92, hep=-0.05", "未达到联合 LOF 阈值"],
      ["MD 轴", "destab=1.158", "tail 接触下降，触发 2M soft evidence"],
      ["当前输出", "2M, conf=0.55", "说明 equilibrium MD face-destab 不是纯 2M 特异"],
    ], pos(70, 160, 790, 294), [140, 210, 440], { rowH: 52, size: 17 });
    addBox(s, pos(910, 170, 240, 240), { fill: "#FFF1F0", line: "#FECACA" });
    addText(s, "错误归因", pos(940, 202, 170, 30), { size: 24, bold: true, color: C.red });
    bulletList(s, [
      "平衡 MD 看不到受力 GOF",
      "contact loss 也可能是 2B release-like",
      "需要 SMD/force-release 轴",
      "2M soft evidence 需加保护条件",
    ], pos(940, 258, 170, 118), { size: 16 });
    addText(s, "优化方向：A1 判断需要 2B_score 与 2M_score 证据竞争。2B_score 不能来自位置记忆，应来自 force-release/SMD、GPIb retained/increased 或独立 GOF 证据。", pos(112, 520, 980, 48), { size: 21, bold: true, color: C.ink });
    foot(s, "No-hotspot Eval: V1316M becomes true 2B → predicted 2M, exposing the need for a directional GOF dynamics axis.");
  }

  // 8
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "系统实现：三层证据进入可解释智能体", "分类器不看 subtype 标签；标签只在预测后用于 Eval。");
    const blocks = [
      ["输入", "• aa_change: L1276P\n• pos=1276, domain=A1\n• ref/alt amino acid"],
      ["特征层", "• AG: RNA/splice delta\n• Boltz: A1-GPIb, heparan\n• MD: contact tail features"],
      ["专家层", "• Structural Expert\n• Transcriptomic Expert\n• Clinical Geneticist"],
      ["融合层", "• interpretable rules\n• 2B/2M directionality\n• uncertain allowed"],
      ["输出", "• subtype: 2M\n• confidence: 0.55\n• reasoning trace"],
    ];
    blocks.forEach((b, i) => {
      const x = 70 + i * 230;
      addBox(s, pos(x, 220, 170, 130), { fill: i === 2 ? "#F7FAFF" : C.light, line: i === 2 ? "#C7D7FE" : C.mid });
      addText(s, b[0], pos(x + 14, 240, 140, 26), { size: 20, bold: true, color: i === 2 ? C.blue : C.ink });
      addText(s, b[1], pos(x + 14, 276, 142, 72), { size: 12, color: C.ink });
      if (i < 4) arrow(s, { x: x + 170, y: 285 }, { x: x + 226, y: 285 });
    });
    addBox(s, pos(118, 452, 1000, 84), { fill: "#F7FFF8", line: "#BBE7C4" });
    addText(s, "无标签泄漏核实", pos(150, 478, 220, 28), { size: 23, bold: true, color: C.green });
    addText(s, "evaluate_vwf_classifier_v2.py 把 true_label/type2_subtype 改成 LEAK_TEST 后，pred_with_md 完全不变；label_leak_test_pass=1。", pos(375, 474, 670, 42), { size: 18 });
    foot(s, "Code refs: scripts/agentic_vwf_classifier.py; scripts/pipeline/evaluate_vwf_classifier_v2.py.");
  }

  // 9
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "Boltz 特征怎么来：输入复合体，输出静态结构置信", "可作为汇报备注：Boltz 不是 MD，它预测一个静态复合体及相关 confidence/score。");
    smallTable(s, ["Boltz 轴", "我们输入什么", "Boltz 输出什么", "我们计算什么"], [
      ["A1-GPIb", "WT/variant A1 + GPIbα 片段", "complex prediction；iPTM/pTM/pLDDT 等", "primary score 与 WT/同轴分布比较 → fb_binding_zscore"],
      ["A1-heparan", "WT/variant A1 + heparan/heparin-like ligand context", "A1 表面结合 context 置信", "同法得到 heparan_zscore；作为 A1 表面/结合面 proxy"],
      ["D'/D3-FVIII", "D'/D3 + FVIII binding context", "复合体 interface confidence", "FVIII 结合轴，用于 2N"],
      ["A3-collagen", "A3 + collagen binding context", "复合体 confidence", "a3_collagen_zscore，用于 A3 型 2M"],
    ], pos(48, 150, 1160, 270), [150, 310, 320, 380], { rowH: 58, size: 13 });
    addBox(s, pos(92, 468, 1030, 88), { fill: "#FFF6F1", line: "#FFD3C2" });
    addText(s, "解释边界", pos(124, 492, 160, 28), { size: 22, bold: true, color: C.orange });
    addText(s, "Boltz 适合问“这个突变是否让某个复合体静态结合/界面置信变差”；不适合直接问“受力后 AIM 会不会打开”或“40-50ns 后接触网络会不会塌陷”。这些属于 MD/SMD 的问题。", pos(300, 486, 760, 46), { size: 17 });
    foot(s, "Sources: Boltz-2 project/technical report; AF3/Boltz confidence metrics used here as proxies, not calibrated binding constants.");
  }

  // 10
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "缺口图谱：哪些能用 AG，哪些不能只靠 AG", "repo 记录：内皮细胞 CL:0000115 有 6 个有效模态；全组织模式可补 5 个表观模态。");
    smallTable(s, ["分型", "现在主要缺什么", "优先补齐工具"], [
      ["Type 1", "定量缺乏：表达/剪接可由 AG 辅助；分泌/清除/ER retention 不能由 AG 直接输出", "AG RNA/splice/regulatory proxy + 结构稳定性/分泌 proxy + 临床 Ag/VWFpp"],
      ["2A", "A2 unfolding、ADAMTS13 可及性、多聚体装配", "AF3/Boltz domain complexes + targeted MD"],
      ["2B", "力依赖 AIM release 还未直接测", "saltbridge retention 仅作 face-retained 线索；SMD/force-MD 补 GOF"],
      ["2M", "GPIb LOF 与 collagen LOF 覆盖不足", "1SQ0 A1-GPIb MD + A3-collagen Boltz/MD"],
      ["2N", "D'/D3-FVIII interface 需要更窄阈值", "Boltz/AF3 complex + PAE/contact calibration"],
      ["Type 3", "null/frameshift/biallelic severe absence 规则", "variant consequence + genotype + AG/NMD"],
    ], pos(50, 150, 1140, 392), [120, 560, 460], { rowH: 58, size: 15 });
    foot(s, "Repo refs: scripts/07e_god_mode_epigenome_crawler.py lines 12-13, 120-122, 494-502; scripts/03_run_alphagenome_inference.py uses CL:0000115 for endothelial RNA/splice.");
  }

  // 11
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "A1 难点：2B 和 2M 都在同一结构域，但方向相反", "2B 是 gain-of-function；2M 是 loss-of-function。去掉位置先验后，这个难点被真实暴露。");
    smallTable(s, ["MD 指标", "直觉解释", "怎么算", "分类含义"], [
      ["saltbridge retention", "AIM-A1 的“锁/桥”是否还保留", "tail window 内关键 AIM-A1 salt bridge/contact 占有率，转 z-score", "高 z：面还在，anti-2M，支持 2B-compatible"],
      ["md_face_destab_score", "A1 结合面是否随时间变不稳", "(0-5ns AIM-A1 masking contacts − 40-50ns tail contacts) / 20", "高：接触面塌陷/不稳，支持 2M LOF"],
      ["tail window", "看生产轨迹最后 40-50ns", "不是每 ns 一个概率分布；是轨迹帧统计后取尾部均值", "减少初始构象/松弛阶段影响"],
    ], pos(48, 148, 1160, 220), [210, 270, 410, 270], { rowH: 60, size: 13 });
    addBox(s, pos(70, 420, 500, 98), { fill: "#F7FFF8", line: "#BBE7C4" });
    addText(s, "新 2M MD 信号", pos(98, 444, 170, 28), { size: 22, bold: true, color: C.green });
    addText(s, "L1276P/R/G1324R/A1377V 从 uncertain → 2M；2M recall +4，但 V1316M 等 2B 会被 face-destab 误伤。", pos(280, 438, 260, 54), { size: 15 });
    addBox(s, pos(640, 420, 500, 98), { fill: "#FFF6F1", line: "#FFD3C2" });
    addText(s, "文献边界", pos(668, 444, 150, 28), { size: 22, bold: true, color: C.orange });
    addText(s, "文献明确 2B/2M 都可涉及 A1-GPIb 轴；2B 常为 A1 GOF/自发 GPIb 结合增强，2M 为 platelet-dependent function 下降。边界不是“域不同”，而是功能方向不同。", pos(820, 430, 290, 66), { size: 14 });
    foot(s, "Sources: ASH/ISTH/NHF/WFH 2021 VWD guideline; Deng Blood 2017 AIM masks A1; PDB 7A6O/Nat Commun 2021 AIM-A1; repo MD report TYPE2M_LOF_MD_FAST_VALIDATION_REPORT.md.");
  }

  // 12
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "MD 与 SMD：先把可用轴跑通", "SMD/force-MD 更接近受力释放机制，但需要慢拉、多重复和严格校准；短期不宜作为第一生产轴。");
    smallTable(s, ["方法", "问的问题", "优点", "瓶颈"], [
      ["Equilibrium MD", "无外力下尾部稳定状态", "可批量；已能提供 2M rescue / anti-2M 线索", "看不到真正剪切/拉伸触发"],
      ["SMD / pull code", "沿反应坐标拉开需要多少力/功", "更接近 AIM release / A2 unfolding", "非平衡；拉太快会失真；需多重复"],
      ["Umbrella/PMF", "自由能曲线", "物理解释更强", "窗口多，成本最高"],
    ], pos(60, 170, 1110, 188), [170, 310, 300, 330], { rowH: 50, size: 15 });
    addBox(s, pos(100, 430, 1000, 92), { fill: "#FFF6F1", line: "#FFD3C2" });
    addText(s, "为什么慢", pos(132, 458, 160, 28), { size: 22, bold: true, color: C.orange });
    addText(s, "50 ns production 已是千万级步数；SMD 还要为每个变体做受力路径、重复轨迹和拉速敏感性分析。GROMACS 文档也强调非平衡拉伸只有在足够慢或多路径统计下才可解释。", pos(300, 450, 740, 48), { size: 17 });
    foot(s, "Sources: GROMACS pull-code docs; VWF A2 force-unfolding MD/AFM literature.");
  }

  // 13
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "三类动力学方法：成本与收益不同", "我们目前做的是第一类 equilibrium MD；SMD/PMF 只适合少量机制验证或高价值 case。");
    const steps = [
      ["Equilibrium MD", "无外力平衡动力学\n50 ns production tail features", "成本：每样本约小时-天级 GPU\n收益：可批量，当前已带来 2B/2M 信号"],
      ["SMD probe", "Steered MD / constant-velocity pulling\n沿反应坐标拉开", "成本：每样本多重复，约数倍-十倍 MD\n收益：可能验证 AIM release / A2 unfolding"],
      ["Umbrella / PMF", "umbrella sampling + potential of mean force\n多窗口自由能曲线", "成本：每样本最高，通常只做少数代表\n收益：机制论文级证据，不适合全量预测"],
    ];
    steps.forEach((st, i) => {
      const x = 84 + i * 380;
      addBox(s, pos(x, 190, 300, 250), { fill: i === 0 ? "#F7FFF8" : i === 1 ? "#FFF6F1" : "#F7FAFF", line: C.mid });
      addText(s, st[0], pos(x + 24, 220, 240, 28), { size: 26, bold: true, color: i === 0 ? C.green : i === 1 ? C.orange : C.blue });
      addText(s, st[1], pos(x + 24, 278, 240, 72), { size: 17 });
      addText(s, st[2], pos(x + 24, 370, 240, 58), { size: 14, color: C.muted });
    });
    addText(s, "汇报口径：SMD/PMF 可作为高价值 case 的机制确认，但短期提升 recall 的性价比仍是 equilibrium MD + score fusion。", pos(112, 520, 990, 36), { size: 21, bold: true });
    foot(s, "Sources: GROMACS pull-code and umbrella sampling docs; force-MD is non-equilibrium and requires replicas/slow pulling for interpretation.");
  }

  // 14
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "预期提升：保持保守，先追求 50-70% recall", "目前不应承诺 80-90%；先把每个分型稳定推入 50-70% 区间更可信。");
    smallTable(s, ["分型", "当前", "补齐机制后目标", "主要抓手"], [
      ["1", "0/115*", "50-60%", "AG RNA/splice/regulatory proxy + secretion/clearance proxy"],
      ["2A", "70/118", "60-70%", "A2/A3/D4 多聚体与 ADAMTS13 轴"],
      ["2B", "15/38", "50-60%", "force-release/SMD + GPIb retained/increased + saltbridge retention"],
      ["2M", "25/53", "50-65%", "A1-GPIb LOF + A3 collagen + MD score fusion"],
      ["2N", "12/16", "60-70%", "D'/D3-FVIII interface tightening"],
      ["3", "未纳入", "60-70% for null-rich", "consequence/genotype/NMD module"],
    ], pos(58, 148, 1120, 362), [120, 160, 220, 620], { rowH: 54, size: 16 });
    addText(s, "* Type 1 当前进入 Eval 但规则未真正覆盖；因此不能当作系统最终能力。", pos(70, 540, 980, 24), { size: 16, color: C.muted });
    foot(s, "Target bands are conservative planning estimates, not validated claims; report them as next-phase milestones.");
  }

  // 15
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "研究定位：不是泛化 pathogenicity，而是 subtype mechanism", "现有 AI 变异工具多回答“是否有害”；我们的问题是“有害通过哪条 VWF 机制表现为哪个 subtype”。");
    smallTable(s, ["相关方向", "已有能力", "留下的空白"], [
      ["AlphaMissense/VEP", "missense pathogenicity prior", "不直接给 VWD subtype 和机制方向"],
      ["AF3/Boltz", "复合体结构与 confidence/affinity proxy", "静态结构不能覆盖动力学/剪切触发"],
      ["Medical LLM agents", "多专家推理框架", "缺少 VWF 特异结构-功能证据链"],
      ["本项目", "结构+AG+MD+规则解释", "需要更完整 Type 1/3 与 2M/2A 机制模块"],
    ], pos(64, 160, 1110, 260), [210, 400, 500], { rowH: 56, size: 16 });
    addBox(s, pos(112, 478, 980, 66), { fill: C.ink, line: C.ink });
    addText(s, "研究问题：能否用结构预测、功能扰动和动力学特征，把 VWD subtype 诊断从表型标签推进为 variant-level 的机制可解释分型？", pos(146, 497, 910, 32), { size: 20, bold: true, color: "#FFFFFF" });
    foot(s, "Sources: AF3 Nature 2024; Boltz-2 technical report/repository; AlphaMissense Science 2023/PubMed.");
  }

  // 16
  {
    const s = deck.slides.add(); s.background.fill = "#FFFFFF";
    title(s, "Take-home", "汇报时建议把“系统还低”说成“机制补全空间明确”。");
    addBox(s, pos(80, 165, 1080, 85), { fill: "#F7FAFF", line: "#C7D7FE" });
    addText(s, "1. 命名可行", pos(112, 190, 180, 30), { size: 24, bold: true, color: C.blue });
    addText(s, "“机制驱动的 VWD 智能体辅助分型/诊断系统”符合项目实质；建议强调 research prototype / decision support。", pos(310, 188, 790, 36), { size: 19 });
    addBox(s, pos(80, 285, 1080, 85), { fill: "#F7FFF8", line: "#BBE7C4" });
    addText(s, "2. 已有价值", pos(112, 310, 180, 30), { size: 24, bold: true, color: C.green });
    addText(s, "删掉 hotspot 后系统更干净；2M MD rescue 更可信，但 2B recall 暴露出缺少 force-release/GOF 动态证据。", pos(310, 308, 790, 36), { size: 19 });
    addBox(s, pos(80, 405, 1080, 85), { fill: "#FFF6F1", line: "#FFD3C2" });
    addText(s, "3. 下一步明确", pos(112, 430, 180, 30), { size: 24, bold: true, color: C.orange });
    addText(s, "补 Type 1/3、重构 A1 score fusion、扩大 2M LOF 与 2A/2N 机制特征，是提升 recall 的主路径。", pos(310, 428, 790, 36), { size: 19 });
    foot(s, "This deck is editable; numbers reflect the latest local Eval artifacts on 2026-06-29.");
  }

  for (const [index, slide] of deck.slides.items.entries()) {
    const stem = `slide-${String(index + 1).padStart(2, "0")}`;
    const png = await deck.export({ slide, format: "png", scale: 1 });
    await writeBlob(path.join(PREVIEW_DIR, `${stem}.png`), png);
    const layout = await slide.export({ format: "layout" });
    await fs.writeFile(path.join(LAYOUT_DIR, `${stem}.layout.json`), await layout.text());
  }
  const montage = await deck.export({ format: "webp", montage: true, scale: 1 });
  await writeBlob(path.join(PREVIEW_DIR, "deck-montage.webp"), montage);
  const pptx = await PresentationFile.exportPptx(deck);
  await pptx.save(FINAL);
  console.log(FINAL);
  console.log(path.join(PREVIEW_DIR, "deck-montage.webp"));
}

main().catch((err) => {
  console.error(err);
  process.exitCode = 1;
});
