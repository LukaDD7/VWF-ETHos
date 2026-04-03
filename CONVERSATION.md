 我看到的项目上下文：

  - CLAUDE.md：/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/CLAUDE.md
    - 说明这是一个基于 AlphaGenome API 的 VWF 变异分析 4 步流程
    - 当前关键脚本是 01_filter_target_vus.py 到 04_analyze_and_visualize.py
    - 记录了核心配置：如 03 用 CL:0000115，04 做 top 20 可视化
    - 还记录了 2026-03-10 的一次完整执行结果：1198 个变异全部跑完并生成图
  - CONVERSATION.md：/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/CONVERSATION.md
    - 这是之前对话/操作日志的转储
    - 重点信息是你最近一直在跑 07e_god_mode_epigenome_crawler_optimized.py
    - 结论里明确写到：
        - 当前 07e 的 有效 success 样本是 145 个
      - 这 145 个样本被认为是 11 模态齐全
      - 一晚运行后没有新增 success，新增的基本都是 VPN 失败记录
      - 失败主因是代理/VPN 流量耗尽，gRPC 连接报错
      - 断点续传/checkpoint/流式写入机制是存在的
  - CHANGELOG*
    - 在这个目录下 没找到 对应文件

  我现在已经有这个项目背景了。接下来你如果要我继续，我可以直接做这几类事之一：

  1. 检查 07e 脚本的失败恢复逻辑是否真的可靠
  2. 重新核对 “145 个 11 模态齐全” 这个结论
  3. 帮你继续分析 07e 结果文件
  4. 看主流程 01-04 代码或输出有没有问题

  你告诉我下一步要查哪一块。

✻ Crunched for 1m 27s

❯ 检查07e，检查主流程01-04的代码逻辑是否有问题，检查当前保存下来的数据有哪些，有什么特点                 

● 我先并行检查代码和结果文件。

  Searched for 5 patterns, read 5 files (ctrl+o to expand)

● 我再核对现有数据文件和关键结果内容。

  Searched for 1 pattern, read 3 files (ctrl+o to expand)

● 我再补几项关键核对：07e结果统计、主流程结果结构、以及07e后续版本是否已经修过部分问题。

  Searched for 2 patterns, read 3 files (ctrl+o to expand)

● 我补一条数量核对。

  Read 1 file (ctrl+o to expand)

● 我看完了。结论先说：

  - 01–04 主流程总体是通的，现有结果也基本自洽
  - 07e 不是“没数据”，而是“有一部分成功，但代码里确实有明显逻辑问题和统计偏差”
  - 当前保存下来的数据很多，最有价值的现成结果是：03/04 主流程全量结果 + 07/07e 的 11 模态整合基础

  ---
  1) 07e 代码逻辑检查结果

  我主要看了：
  - scripts/07e_god_mode_epigenome_crawler_optimized.py
  - scripts/07e_god_mode_epigenome_crawler_optimized_v2.py

  07e 的主要问题

  问题 A：success_count = len(completed_variants) 统计是错的

  位置：
  - scripts/07e_god_mode_epigenome_crawler_optimized.py:435
  - scripts/07e_god_mode_epigenome_crawler_optimized_v2.py:603

  这里把“已完成变异数”直接当成“成功数”。

  但 completed_variants 里其实同时包含：
  - Success
  - Failed
  - Pending

  来源在：
  - scripts/07e_god_mode_epigenome_crawler_optimized.py:174-177
  - scripts/07e_god_mode_epigenome_crawler_optimized_v2.py:283-286

  所以这会导致：
  - 日志里的成功计数不可信
  - 进度显示会混淆“已处理”和“成功”
  - 你之前看到的一些“成功率/剩余量”描述可能被误导

  ---
  问题 B：断点续传把 Failed 也当成“已完成并永久跳过”

  位置：
  - scripts/07e_god_mode_epigenome_crawler_optimized.py:174-177
  - scripts/07e_god_mode_epigenome_crawler_optimized.py:445-449

  当前逻辑是：
  - 只要 CSV 里这个位点出现过，且状态属于 Success/Failed/Pending
  - 下次就直接跳过

  这意味着：
  - 网络失败的样本不会自动重试
  - 真正的坏样本和 VPN 失败样本被混在一起永久跳过
  - 这和“断点续传自动重跑失败项”的直觉并不一致

  这也是 07e 最核心的问题之一。

  ---
  问题 C：变异唯一键只用 chrom:pos，没有带 ref/alt

  位置：
  - scripts/07e_god_mode_epigenome_crawler_optimized.py:176
  - scripts/07e_god_mode_epigenome_crawler_optimized.py:442
  - v2 同样如此

  现在 key 是：
  - chr12:6019618

  但同一位置可能有多个 ALT，甚至不同 REF/ALT 组合。
  你当前数据里就能看到同位点多变异，例如 checkpoint 里出现：
  - chr12:6019618 A>T
  - chr12:6019618 A>G

  如果只用坐标当 key，会造成：
  - 一个变异跑过后，同位点另一个变异可能被错误跳过
  - 结果数和任务数的映射失真
  - “11 模态齐全有多少” 这类统计容易被污染

  这个问题很严重。

  ---
  问题 D：get_max_delta_with_idx() 取索引的实现有 bug，解释了你看到的大量 index ... out of bounds

  位置：
  - scripts/07e_god_mode_epigenome_crawler_optimized.py:291-310
  - v2 同样如此

  这里对多维数组用了：

  max_idx = int(np.argmax(delta))
  max_val = float(delta.flat[max_idx] if delta.ndim == 1 else delta[..., max_idx].max())

  对于二维/多维数组：
  - np.argmax(delta) 返回的是 flatten 后的一维索引
  - 但后面又把这个 flatten 索引当成了最后一维 track 索引去切 delta[..., max_idx]

  这会直接导致你日志里那种错误：

  - index 695216911 is out of bounds for axis 1 with size 667
  - index 627461717 is out of bounds for axis 1 with size 667

  也就是说：
  - 很多 07e 失败不是 API 真失败
  - 而是脚本自己在后处理 delta 时索引错了

  这是 07e 里最明确的代码 bug。

  ---
  问题 E：1Mb 区间构造没有像 02 那样做染色体边界保护

  位置：
  - scripts/07e_god_mode_epigenome_crawler_optimized.py:356-357

  interval_start = max(1, pos - half_len)
  interval_end = pos + half_len

  这里只保护了 start，没有保护 end。
  虽然 VWF 在 chr12 中段，当前这批位点大概率不会越界，但逻辑上不完整。

  ---
  问题 F：失败类型没有区分“网络失败”和“数据本身非法”

  你 checkpoint 里已经能看到两类失败：

  1. 网络失败
    - failed to connect to all addresses
    - recvmsg: Connection reset by peer
  2. 输入数据/脚本失败
    - Invalid reference bases: "c". Must only contain "ACGTN".
    - index ... out of bounds for axis 1 with size 667

  现在它们都被统一记成 Failed，然后都在断点续传时跳过。
  这会让后续恢复很难。

  ---
  07e 当前结果到底是什么状态

  我核对了文件：

  - results/07e_GodMode_Epigenome_Peaks.csv：506 行
    - 1 行表头
    - 505 条记录
  - 其中：
    - 145 条 Success
    - 193 条 Failed
  - 但 CSV 总记录数 505，不等于 145 + 193 = 338

  这说明：
  - 文件里还存在一批 非 ,Success / ,Failed 结尾的记录
  - 很可能是格式受损、状态异常、或历史版本残留行
  - 再结合 checkpoint 只有 370 行，CSV 与 checkpoint 也不同步

  所以 07e 当前结果文件不是完全干净的审计真源。
  更可信的是：
  - 成功样本数 145 这个结论目前能确认
  - 但“失败多少”“处理了多少”“还能续跑多少”需要重新清洗后才能精确说

  ---
  2) 主流程 01–04 代码逻辑检查结果

  我看了：
  - scripts/01_filter_target_vus.py
  - scripts/02_preprocess_and_liftover.py
  - scripts/03_run_alphagenome_inference.py
  - scripts/04_analyze_and_visualize.py

  01_filter_target_vus.py

  整体逻辑是合理的。

  优点：
  - 会检查输入文件存在
  - 会检查关键列存在
  - 同时保留 Uncertain_significance + Pathogenic + Likely_pathogenic
  - 把致病变异保留作阳性对照，这个思路没问题

  注意点：
  - 文件头注释写的是“过滤目标 VUS”，但实际代码保留了致病/可能致病
  也就是 注释比代码更窄，代码行为更广
  - matches_consequences() 用的是子串匹配，能工作，但比较宽松
  不过在你这类 consequence 复合字段里，实际上是可接受的

  结论：
  - 没发现致命逻辑 bug
  - 只是命名/注释与实际策略不完全一致

  ---
  02_preprocess_and_liftover.py

  整体也没大问题。

  优点：
  - 使用本地 chain file，流程稳
  - 对 liftover 失败做了严格剔除
  - 构建 1Mb interval
  - 对染色体边界做了修正

  注意点：
  - liftover_hg19_to_hg38() 直接取 result[0]，默认取第一条映射
  对大多数 VWF 位点没问题，但理论上多重映射没被特别处理
  - 输出列有点保守，只保留了一部分原始字段
  如果后续想追更多注释，需要回到原表

  结论：
  - 逻辑基本正确
  - 没看到会系统性污染结果的大 bug

  ---
  03_run_alphagenome_inference.py

  这个脚本总体也比较稳，且结果和现有输出是对得上的。

  优点：
  - 有 checkpoint：03_inference_checkpoint.jsonl
  - 成功/失败都会记录
  - 结束时会合并输出 CSV + PKL
  - 优先处理 pathogenic，当作阳性对照，这个合理

  注意点：

  问题 1：API key 明文写在代码里

  位置：
  - scripts/03_run_alphagenome_inference.py:46

  这不是逻辑 bug，但属于安全/维护问题。建议改成环境变量。

  问题 2：checkpoint 恢复后的 raw_outputs 丢失

  位置：
  - scripts/03_run_alphagenome_inference.py:287-324

  从 checkpoint 恢复时，只恢复摘要字段，raw_outputs 没法恢复。
  这意味着如果中途中断后再导出最终 PKL，checkpoint 恢复出来的历史结果不会带完整 raw object。

  不过你当前的 03_inference_results.pkl 已经是完整成品，体积 57G，说明那次全量成功运行时没出这个问题。

  问题 3：如果输入里有非法碱基，没有预清洗

  但主流程现有结果 1198 条都是成功的，所以这批 03 输入质量应该是好的。

  结论：
  - 主流程脚本 03 总体可用
  - 真正的大成品文件也说明这套流程跑通了
  - 最主要的问题是 checkpoint 恢复后不适合再做需要 raw_outputs 的完整 PKL 重建

  ---
  04_analyze_and_visualize.py

  总体逻辑可以，且现有输出与之匹配。

  优点：
  - 从 03_inference_results.pkl 读取 raw outputs
  - 按 delta_score 排前 TOP_N_VARIANTS
  - 画 RNA / splice 双模态 4 面板图
  - 输出 04_analysis_summary.csv

  注意点：

  问题 1：强依赖 03_inference_results.pkl 中存在 raw_outputs

  位置：
  - scripts/04_analyze_and_visualize.py:77
  - scripts/04_analyze_and_visualize.py:186

  如果 03 是靠 checkpoint 恢复重建出来的 PKL，这里可能会失效。
  但你当前这个 57G PKL 是完整的，所以现在没事。

  问题 2：CENTER_IDX = 524288 是硬编码

  位置：
  - scripts/04_analyze_and_visualize.py:92

  这依赖于输入永远是 1Mb 且变异居中。
  在你当前流程里是成立的。

  结论：
  - 04 逻辑没明显问题
  - 现有的 04_analysis_summary.csv 看起来也正常

  ---
  3) 当前保存下来的数据有哪些

  data/

  主要是主流程输入中间件：

  - data/01_filtered_targets.csv — 1.8M
  - data/02_ready_for_inference.csv — 1.4M
  - data/Clinvar_HGMD_merge_Hg37_fixed.vcf
  - data/Clinvar_HGMD_merge+分析20260312_删除32个不能注释.xlsx
  - data/hg19ToHg38.over.chain.gz

  说明：
  - 主流程 01、02 的中间产物都在
  - liftover 依赖文件也在

  ---
  results/

  这是最重要的目录，体量很大：609G

  主流程 01–04 成果

  - results/03_inference_results.csv — 1.4M
  - results/03_inference_results.pkl — 57G
  - results/04_analysis_summary.csv — 1.4M

  数量上：
  - 03_inference_results.csv：1199 行 = 1198 条结果 + 表头
  - 04_analysis_summary.csv：1199 行 = 1198 条成功排序结果 + 表头

  这和 CLAUDE 里记的“1198/1198 成功”一致。

  07 / 11模态相关大文件

  - results/07_VCF_AlphaGenome_Results.csv — 313K
  - results/07_VCF_AlphaGenome_Results.pkl — 185G
  - results/07_VCF_AlphaGenome_Results_Backup.csv — 313K
  - results/07_VCF_AlphaGenome_Results_Backup.pkl — 185G
  - results/07_VCF_AlphaGenome_Results_11_Modalities_Final.pkl — 185G

  说明：
  - 你已经有一套非常大的 07 系列整合结果
  - CSV 轻量版可以做任务索引，PKL 是全量对象仓库
  - 11_Modalities_Final.pkl 很可能是你后续核心资产之一

  07e God Mode 相关

  - results/07e_task_list.json — 233K
  - results/07e_GodMode_Epigenome_Peaks.csv — 43K
  - results/07e_GodMode_Epigenome_Peaks_backup_20260324_155250.csv — 571K
  - results/07e_godmode_checkpoint.jsonl — 188K

  说明：
  - 07e 当前成品是 轻量补充表
  - 但是现有 CSV/JSONL 状态不够干净，需要重新审计

  最终交付类表格

  - Final_VWF_Lossless_Merge.xlsx
  - Final_VWF_Target_List_with_AlphaGenome.xlsx
  - Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx

  ---
  4) 当前数据的主要特点

  主流程 03/04 的特点

  从 results/03_inference_results.csv 和 results/04_analysis_summary.csv 看：

  - 1198 个变异都有结果
  - 全部集中在 chr12 / VWF
  - 排名前列的变异，delta_score 基本由 RNA delta 主导
  - splice_site_delta_max 常常接近 1，但不是排序主导项
  - top hits 里富集：
    - frameshift
    - stop_gained
    - disruptive_inframe_deletion
    - splice-related variants

  这和生物学直觉是吻合的。

  ---
  07e 成功数据的特点

  从 results/07e_GodMode_Epigenome_Peaks.csv 前 20 行看：

  - 当前成功样本都在 chr12
  - 成功记录里 5 个补充模态都有值：
    - ATAC
    - Histone
    - CAGE
    - PROCAP
    - Contact
  - Contact_Max_Delta 前面大量是 0.0625
  - RNA_Delta_HUVEC 常见值：
    - 2.5
    - 4.0
  - Splice_Delta_HUVEC 常见值很小：
    - 0.01171875
    - 0.015625
    - 0.03125
  - Histone / CAGE 常是比较大的离散值：
    - 512, 640, 768, 896 等

  这提示：
  - 07e 成功样本不是空值，它们确实有稳定补充模态信号
  - 但这些值里有明显的离散/量化痕迹，后续分析时要注意不要把它们当成连续精细量纲

  ---
  07e 失败数据的特点

  checkpoint 里可见三类失败：

  1. VPN / 代理失败
    - StatusCode.UNAVAILABLE
    - recvmsg: Connection reset by peer
  2. 输入碱基格式问题
    - 例如 ref = "c" 小写
    - AlphaGenome 要求 ACGTN
  3. 代码后处理索引 bug
    - index ... out of bounds for axis 1 with size 667

  也就是说：
  - 07e 当前失败并不都代表样本无效
  - 很大一部分只是：
    - 网络失败
    - 代码 bug
    - 输入规范化不足

  ---
  5) 最重要的判断

  可以信的

  - 主流程 01–04 的现有成品基本可信
  - 03_inference_results.csv / 04_analysis_summary.csv 是可直接分析的
  - 07e 已成功拿到 145 个补充 5 模态样本

  不能直接信的

  - 07e 里的失败数、处理进度、自动续跑状态
  - “Failed 就代表样本本身不行” 这个结论
  - 07e 现有断点续传逻辑

  ---
  6) 我建议的下一步

  最值得做的是这 3 件事，按优先级排：

  方案 1：先修 07e 逻辑

  重点修：
  1. variant key 改成 chrom:pos:ref:alt
  2. Failed 不要自动视为永久完成
  3. 修 get_max_delta_with_idx() 的索引 bug
  4. 输入 ref/alt 统一大写
  5. 分开记录 network_failed / invalid_variant / extract_bug

  方案 2：清洗当前 07e 结果

  我可以帮你直接产出一份审计表，分成：
  - 真成功
  - 网络失败待重跑
  - 输入非法待修正
  - 脚本 bug 导致待重跑
  - 重复位点/重复变异

  方案 3：核对“11 模态齐全”的严格定义

  因为 07e 现在 key 有问题，我建议重新严格核对一次：
  - 按 chr:pos:ref:alt
  - 和 07_VCF_AlphaGenome_Results_Backup.csv 做精确匹配
  - 算出真实的 11 模态齐全数量

  ---
  如果你愿意，我下一步可以直接做这两件事之一：

  1. 修 07e 脚本
  2. 先把当前 07e 数据做一次彻底审计统计

  我建议先做 第 2 个，因为先把脏数据分层，你再决定要不要重跑。


  已在 /media/luzhenyang/project/alphagenome/alphag
  enome/VWF_ErTongYiyuan 目录。

  该目录包含以下内容：

  主要数据文件：
  - 07_VCF_AlphaGenome_Results_11_Modalities_Final.
  pkl (180GB) - 完整的11模态AlphaGenome分析结果
  - VWF_Type2_variants.csv - Type 2 VWF变异列表
  - VWF_Type2_AF3_Reference_Table.xlsx/csv -
  AlphaFold3参考表

  分析脚本：
  - vwf_type2_analysis.py - Type 2变异分析主脚本

  子目录：
  - Proteo-Structure-Pipeline/ -
  AlphaFold3蛋白结构预测流程
  - scripts/ - 各类分析脚本
  - figures/ - 可视化输出
  - results/ - 分析结果
  - structures/ - 结构预测结果
  - data/ - 数据文件
  - logs/ - 日志文件

  文档：
  - README.md - 项目说明
  - Type-2-Analysis-README.md - Type 2分析详细说明
  - 研究思路.md - 研究思路文档

!cd /media/luzhenyang/project/alphagenome/alphagen 
  ⎿  (Bash completed with no output)               

❯ 去看看type-2型的AF3结果有多少了,然后简单统计下具 
体分型的分布/media/luzhenyang/project/alphagenome/ 
alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipe 
line/output/af3_batches_type2/AF3_Results          
  Read 2 files (ctrl+o to expand)

● Bash(cd /media/luzhenyang/project/alphagenome/alp
      hagenome/VWF_ErTongYiyuan/Proteo-Structure-Pi
      peline/output/af3_batches_type2/AF3_Results
  ⎿  64& ls extracted/ | wc -l &&…)
     fold_vwf_a1437t
     fold_vwf_a1461v
     … +18 lines (ctrl+o to expand)

● Bash(cd /media/luzhenyang/project/alphagenome/alp
      hagenome/VWF_ErTongYiyuan/Proteo-Structure-Pi
      peline/output/af3_batches_type2/AF3_Results
  ⎿  总用量 248 extracted/ | head…)
     drwxrwxr-x. 66 luzhenyang luzhenyang 4096 3月  22 15:04 .
     drwxrwxr-x.  5 luzhenyang luzhenyang 4096 3月  29 11:36 ..
     … +26 lines (ctrl+o to expand)

● Bash(python3 << 'EOF'
      import pandas as pd…)
  ⎿  ============================================================
     Type-2 VWF变异 AF3结构预测结果统计
     ============================================================
     … +48 lines (ctrl+o to expand)

● Type-2型VWF变异AF3结构预测结果统计

  完成情况

  - 已提取结构数: 64个目录 (包含59个变异 + WT对照)
  - 已完成分析: 59个Type-2变异
  - 平均pLDDT: 74.13 ± 0.18 (高质量预测)

  ---
  1. Type-2分型分布 (59个变异)

  ┌────────┬──────┬───────┬────────────────────────┐
  │  分型  │ 数量 │ 占比  │          说明          │
  ├────────┼──────┼───────┼────────────────────────┤
  │ type2A │ 36   │ 61.0% │ 最常见，功能缺陷型     │
  ├────────┼──────┼───────┼────────────────────────┤
  │ type2N │ 13   │ 22.0% │ 抗原缺乏型             │
  ├────────┼──────┼───────┼────────────────────────┤
  │ type2M │ 8    │ 13.6% │ 功能缺陷型(多发于A3域) │
  ├────────┼──────┼───────┼────────────────────────┤
  │ type2B │ 2    │ 3.4%  │ 功能增强型(A1域)       │
  └────────┴──────┴───────┴────────────────────────┘

  ---
  2. ACMG致病性分类

  ┌──────────────────────────────┬──────┬───────┐
  │             分类             │ 数量 │ 占比  │
  ├──────────────────────────────┼──────┼───────┤
  │ Uncertain significance (VUS) │ 23   │ 39.0% │
  ├──────────────────────────────┼──────┼───────┤
  │ Pathogenic                   │ 18   │ 30.5% │
  ├──────────────────────────────┼──────┼───────┤
  │ Likely pathogenic            │ 18   │ 30.5% │
  └──────────────────────────────┴──────┴───────┘

  ---
  3. VWF功能域分布

  ┌───────────────────────┬──────┬───────┬──────────────────────────────────┐
  │        功能域         │ 数量 │ 占比  │             主要分型             │
  ├───────────────────────┼──────┼───────┼──────────────────────────────────┤
  │ A1 (胶原结合，type    │ 19   │ 32.2% │ type2N (13), type2A (6)          │
  │ 2B)                   │      │       │                                  │
  ├───────────────────────┼──────┼───────┼──────────────────────────────────┤
  │ A3 (胶原结合，type    │ 15   │ 25.4% │ type2A (12), type2M (3)          │
  │ 2M)                   │      │       │                                  │
  ├───────────────────────┼──────┼───────┼──────────────────────────────────┤
  │ D1-D2 (信号肽+前肽)   │ 7    │ 11.9% │ type2A (7)                       │
  ├───────────────────────┼──────┼───────┼──────────────────────────────────┤
  │ A2 (ADAMTS13切割)     │ 7    │ 11.9% │ type2M (4), type2B (1), type2A   │
  │                       │      │       │ (2)                              │
  ├───────────────────────┼──────┼───────┼──────────────────────────────────┤
  │ D4 (多聚化)           │ 5    │ 8.5%  │ type2M (2), type2A (3)           │
  ├───────────────────────┼──────┼───────┼──────────────────────────────────┤
  │ D'D3 (GP1b结合)       │ 4    │ 6.8%  │ type2A (4)                       │
  ├───────────────────────┼──────┼───────┼──────────────────────────────────┤
  │ CT (C端)              │ 2    │ 3.4%  │ type2A (2)                       │
  └───────────────────────┴──────┴───────┴──────────────────────────────────┘

  ---
  4. 分型 × ACMG交叉统计

  ┌────────┬───────────────────┬────────────┬─────┬──────┐
  │  分型  │ Likely Pathogenic │ Pathogenic │ VUS │ 总计 │
  ├────────┼───────────────────┼────────────┼─────┼──────┤
  │ type2A │ 13                │ 8          │ 15  │ 36   │
  ├────────┼───────────────────┼────────────┼─────┼──────┤
  │ type2N │ 3                 │ 7          │ 3   │ 13   │
  ├────────┼───────────────────┼────────────┼─────┼──────┤
  │ type2M │ 1                 │ 3          │ 4   │ 8    │
  ├────────┼───────────────────┼────────────┼─────┼──────┤
  │ type2B │ 1                 │ 0          │ 1   │ 2    │
  └────────┴───────────────────┴────────────┴─────┴──────┘

  关键观察:
  - type2N: 致病性变异比例最高 (10/13 = 77%)
  - type2A: VUS比例最高 (15/36 = 42%)
  - type2M: VUS比例较高 (4/8 = 50%)
  - type2B: 样本最少，致病性和VUS各半

  ---
  5. 各分型结构质量

  所有分型pLDDT均在74左右，质量稳定:
  - type2A: 74.14 (73.69-74.63)
  - type2N: 74.14 (73.89-74.32)
  - type2M: 74.09 (73.81-74.36)
  - type2B: 74.09 (74.02-74.16)

  高置信度原子比例均~72%，预测质量优良。