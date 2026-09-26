# 实现契约（格式 / 解析 / 测试约束）

> **本文件是工程约束，不是生物学判据。** 生物证据规则见 `evidence-standard.md` / `annotation_quality.md`；
> 这里的内容回答"代码必须怎样表现、由哪组测试锁住"。
> 拆出本文件的原因：把"selector 归一化后为空则加载失败"这类**实现约束**写在生物学文件里，
> 容易被后来者误读成"领域规则"。

## 1. 两条验证路径必须等价

案例记录有两条校验路径：`schemas/case.schema.json` + `jsonschema`，以及无依赖时的显式兜底
`_case_errors_without_jsonschema()`。**两者必须对同一份 schema 给出一致裁定**，由固定数据集锁住：

| 测试 | 锁住什么 |
|---|---|
| `tests/test_evidence_contracts.py::SchemaValidatorEquivalenceTests` | 合法 / 缺字段 / 类型错 / 非法枚举 / 嵌套对象错 / 跨字段冲突 → 两路径裁定相同 |
| `tests/test_pr1_review_regressions.py::SchemaKeywordCoverageTests` | schema 的每个顶层关键字至少有一条"兜底也必须拒绝"的样本 |
| `tests/test_pr1_review_regressions.py::...fallback_verdicts_are_fixed` | 无 `jsonschema` 环境下的裁定快照 |

**改动 schema 时必须同时跑这两组测试**：否则两条路径会再次静默偏离（历史上已经偏离过两次：
显式 `decision: null`；`case_id` 的 `minLength: 1` 与三个可选数组）。
CI 固定安装 `jsonschema`，因此等价性测试在 CI 中真正执行而不是跳过。

**跨字段业务规则不在 schema 里表达**（如案例级 `RESOLVED` 与某异常 `UNRESOLVED` 并存）；
若将来要加，必须**同时**加在两条路径上（见 §7 技术债）。

## 2. `/transl_except` 与例外 registry 的加载约束（fail-closed）

生物学含义（何时算已验证）见 `annotation_quality.md` §9；**格式非法即受控失败**的清单如下：

1. **selector 归一化后为空或显式为 `null` → 加载失败**：`gene` 为 `""`/`"?"`/纯标点/`"(CUN)"`，
   `codon`/`amino_acid` 为空字符串，或键存在但值为 `null`。
   理由：空 selector 会与"缺 `/gene`+`/product` 的 CDS"（canonical `""`）**精确匹配**，
   等于把已审计例外挂到"没有任何基因身份"的记录上。
2. **记录类型按"键是否存在"判定**（而非真值）：`amino_acid` 键一旦出现就是 `transl_except` 记录，
   `""` 或 `null` 均为格式错误，**不得降格**为 start 记录。
3. **`transl_except` 的 `codon`** 必须是三个 IUPAC 碱基；`amino_acid` 必须是合法例外 token（不得为 `TERM`）。
4. **start selector 不得为空**：`--tolerate-start` 的基因名或 start 记录的 `gene` 归一化后为空
   （`":CGA"`、`gene:"?"`/`""`/纯标点）→ 参数/加载受控失败。
5. **统一密码子规则**（start 记录、`transl_except` 记录、`--tolerate-start` 共用同一函数）：
   先 `.strip().upper()`，再要求**恰好三个 IUPAC 碱基**（`[ACGTURYKMSWBDHVN]{3}`）。
   `""` / `"   "` / `"CG"` / `"XXXX"` / `"C1A"` 在加载/参数阶段受控失败；`"cga"` 与 `" CGA "` 归一化为 `"CGA"`。
6. **start 记录按 `(gene, codon)` 唯一**：重复（含仅 `taxon` 不同）→ **加载失败**，不允许后一条静默覆盖；
   若确需多类群共存，必须先扩展 start registry 结构并定义选择规则。
7. **`transl_except` 记录同位点、不同 `taxon`/`transl_table` 可共存**（record list + 显式选择）；
   仅当 selector 全等（gene/codon/amino_acid/位点/taxon/transl_table）才判重复并拒绝。
8. **位点必须绑定**：记录需含 `codon_index` 或 `pos`，或显式 `scope="gene_wide"`；
   `transl_table`/`taxon`/`source`/`rationale` 必须为正确类型且非空。
9. **未登记/未验证只维持 REVIEW**，不升级为已验证：`EXCEPTION_NOT_REGISTERED` / `EXCEPTION_RECORD_INCOMPLETE` /
   `EXCEPTION_TAXON_MISMATCH` / `EXCEPTION_TAXON_UNVERIFIED`；
   而**格式非法（本节第 1–8 项）导致加载失败并退出 1**，不进入 REVIEW，也不做部分生效。

## 3. 退出码语义（不得误读）

| 工具 | 码 | 含义 | **不是** |
|---|---|---|---|
| `annot_check.py` | `0` / `1` / `2` | 无发现 / 有错误（含参数错误）/ 仅待核查 | `2` **不是**通过；`case-validate` 的通过也不是科学验收 |
| `cox1_id.py` | `0` / `1` / `2` / `3` | 得判读 / `insufficient` 或 `no_match` / 拒绝执行 / 网络或结果格式故障 | `3` **不是**"没找到相似序列" |
| `blast_genes.py` | `0` / `1` | 全部基因唯一命中 / 任一基因无唯一命中或命中歧义（**不写部分结果**） | `1` 不是"基因真的缺失" |
| `mitos2_to_genbank.py` | `0` / 非 0 | 写出 `<out.gb>` / 缺文件、坐标越界、`result.fas` 与坐标不一致、未知 feature 类型 | 非 0 时**不得**留下半成品 |
| `depth_analysis.py` | `0` / `1` / `2` | 正常 / 输入或 BAM 错误 / 存在低覆盖区 | `2` 不是"组装错误" |
| `circularize.py` | `0` / `2` / `1` | 接受候选 / REVIEW（含待回贴）/ 失败或接缝证据不足 | 任何一个码都不能独自证明物理闭环 |

后台任务"子进程正常退出"只表示命令执行成功；科学验收由证据标准另行判定。

## 4. 写入安全（工程要求，与生物学无关）

- 候选修复写入**独立目录**，绝不覆盖原始输入；
- 结构化案例的写入是**原子**的（临时文件 + `os.replace`），失败不留 `.tmp`、不留半成品 `case.json`；
- `case-anomaly` 的重复 `id` **fail closed**（文件字节不变，`--update` 才替换）；
- `--event-action` 必须命中 `events.jsonl` 中真实存在的事件（避免悬空关联）；
- `mitos2_to_genbank.py` 的拓扑参数只写声明，**MITOS2 的 circular 模式不是物理环化证据**。

## 5. 测试与规则的对应关系

失败复现测试**永久保留**；每条硬规则都要有对应的锁：

| 规则 | 锁定测试 |
|---|---|
| GenBank 语义（链、跨原点、`transl_except`、partial） | `tests/test_annot_semantics.py`、`tests/test_annot_check_fixes.py` |
| COX1 结构化 BLAST（XML2/旧 XML/半可读 HSP=格式故障） | `tests/test_annot_cox1_edges.py`、`tests/test_pr1_review_regressions.py` |
| 两条校验路径等价 | `tests/test_evidence_contracts.py`、`tests/test_pr1_review_regressions.py` |
| 工具链桥（MITOS2→GenBank、`run_mitos2` outdir、`check_env` 提示） | `tests/test_toolchain_fixes.py` |
| 案例 CLI（`case-anomaly`、假设编号） | `tests/test_case_anomaly_cli.py` |
| 公共参考登记（accession.version、等级→用途、漂移、截断） | `tests/test_reference_registry.py` |
| 案例类型与 lesson 范围/领域（复审核 R-1…R-8） | `tests/test_case_type_and_lesson_scope.py`、`tests/test_review_round_lesson_domain_and_report.py` |
| 注释策略与类群例外 | `tests/test_annotation_policy.py` |
| 数据与路径保护（未发表数据不进公共库） | `.gitignore` 规则 + `git check-ignore`/`git add -n` 人工核验（见 `REAL-DATA` 审计记录） |

新增/修改行为时的顺序：**先写失败复现 → 最小修复 → 补反例 → 跑完整测试与 CI**，CI 绿灯是必要条件。

## 6. 文档分层约定（避免把实现约束写成领域规则）

| 层 | 归属 | 例子 |
|---|---|---|
| (a) 一般生物学预期 | `annotation_quality.md` / `standard_gene_order.md`（须附来源与覆盖类群） | 13 CDS / 22 tRNA 是**预期不是定律** |
| (b) NCBI 提交审查要求 | `evidence-standard.md` §6（须引用官方页面） | 提交时需向策展人说明差异 |
| (c) 本工具工程阈值 | 相关文件显式标注 | `>8bp` 重叠、tRNA `60–75bp`、`9+/4−`、`max_span_ratio=3.0` |
| (d) 实现约束 | **本文件** | selector 归一化、原子写入、退出码、schema 等价性 |

(c)/(d) 层**不得**写成"NCBI 要求…"或领域公理。

## 7. 已知技术债

1. **两套案例校验实现**（§1）：保持等价是硬要求，但仍是重复实现。
2. **跨字段矛盾不校验**（如案例级 `RESOLVED` 与异常 `UNRESOLVED` 并存）：schema 不表达，
   两套实现都不拒绝；若要加必须同时加在两边。
3. **BLAST 坐标边界**：`cox1_id.py` 的 HSP 边界情形仍有未覆盖组合（详见 issue 记录）。
4. **`case_type` 是可选字段**：`learning-policy.md` 把记录分类列为必需，但 schema 未列入 `required`。
   若改成强制项，已有/外部归档案例（可能缺该字段）会一次性变成 `INVALID`，
   因此当前策略是**可选 + 缺失时按 `abnormal_case` 保守处理 + 新案例一律写入**；
   要真正强制，必须同时提供迁移方案（回填 `case_type`）与两套校验路径的同步修改。
5. **lesson 的领域/范围不是 schema 强制项**：写入端（`propose-lesson`）是 fail-closed 的，
   验收端（`validate_public_lesson`）对缺失值做保守补齐后才接受；两边规则不同步就会再造出偏差，
   因此改任意一边必须同时跑 `tests/test_review_round_lesson_domain_and_report.py`。

6. **报告按状态机械分类**：`case_report()` 将 RESOLVED/NO_CHANGE 放入“有证据支持”，
   UNRESOLVED 放入“无证据支持”。这不符合证据语义；当前只能把输出当草稿，按
   [conclusion-report.md](conclusion-report.md) 复核。下一步修复应以部分支持但未解决、
   证据冲突、无证据的 NO_CHANGE 等案例验证，不能用改标题掩盖证据欠缺。
7. **科学报告字段未完整承载**：schema/CLI 未完整校验支持/反证、未测项、范围、局限、
   结论类型与复核信息。`--event-action` 是可选的，只验证被传入的 action 存在；
   action 同名时无法唯一定位某次执行。需由事件与报告补足，不能称“格式通过即完整可追溯”。
   `case-anomaly --update` 会替换该条记录；手工附加字段不能依赖它保留。
8. **类群配置与经验独立性**：`--taxon` 不自动切换昆虫阈值；lesson 的“不同 case_id
   且不同类群”是当前推广门槛的工程计数，不等于科学独立性。文档解释层不得冒称代码已修复。
9. **覆盖提示过强**：`depth_analysis.py` 会将高 soft-clip 比例按区域标签提示为“正常长度异质性”
   或“缺失序列”。这些是未验证的解释，不能作为诊断事实；当前按证据标准重新判读。

## 8. 注释解析细节（仅调试或维护时读取）

以下保留原注释手册中的解析/格式约束；科学接受条件仍见
[annotation_quality.md](annotation_quality.md) §2/§9。下面的 JSON 仅演示格式，含占位来源，
不是可直接使用的已审计生物学例外。

   - `/transl_except` **分三层，存在不等于豁免**：
     **① 语法层**：代码按 CDS 转录顺序建立每个密码子的**精确基因组位置集合**；整个 qualifier 必须被
     **完整消费**（一个或多个 `(pos:...,aa:...)`），语法不完整就记 `TRANSL_EXCEPT_UNPARSED` 并**整条作废，
     不沿用合法前缀**（尾随逗号、尾随损坏文本、缺 `pos:` token、模糊位置、括号不配对都属此类）。
     `pos` 支持 `a..b`、`complement(a..b)`、`join(..)`/`order(..)`（密码子可跨 `join()` 边界）；
     同一 `pos` 内**不同链方向的片段混用直接拒绝**（整段位置不得由 OR 合并成单链）。
     `aa` 须是合法三字母代码且**不是 `TERM`**；`pos` 方向须与 CDS 一致（**负链 CDS 必须 `pos:complement(a..b)`**，
     正链不得写 `complement`）；位置集合必须**恰好**等于某个真实内部终止密码子。通过只记
     `TRANSL_EXCEPT_MATCHED`（位置/读框事实），**MATCHED 本身不等于已接受**。
     **② 生物学层**：语法匹配后还必须有审计证据才能接受。证据来自 `--exception-registry` 的
     `transl_except` 条目，键为 `gene + codon + amino_acid`，且必须满足三个 fail-closed 约束：
     （a）**样本绑定**：必须显式给出 `--taxon` 并与记录 `taxon` 一致 —— 未提供 `--taxon` 时无法把记录绑定到本样本，
     **不得**升为已验证；
     （b）**位点绑定**：记录必须说明它覆盖**哪个** stop —— 给 `codon_index`（CDS 转录顺序上的 1-based 密码子序号）
     或 `pos`（该密码子的基因组位置，如 `10..12`）之一；只有**显式** `scope="gene_wide"` 的记录才允许复用多个位点。
     不绑定位点的记录在**加载时**受控失败（否则一条记录会静默变成基因范围的重编码规则）；
     （c）**字段完整且类型正确**：`transl_table` 须为正整数且等于 `--table`；`source`、`rationale` 须为非空字符串；
     数组/对象/整数等错误 JSON 类型在**加载时**受控失败（`str(value)` 非空不等于已审计）。
     三者均满足才记 `TRANSL_EXCEPT_VALIDATED`（并输出 `scope=`）并扣除该项。
     **不引入全局密码子→氨基酸重编码表**：密码子含义随类群与密码表变化，把“合法 INSDC token”当成
     “已验证例外”会把类群特异的例外伪装成普遍规律。因此 `TAA -> Gln` 这类任意声明即使位置匹配，
     也只能得到 `TRANSL_EXCEPT_DECLARED_UNVERIFIED`（REVIEW），**内部终止仍保留 ERROR**。
     **③ 结果层**：因此**一条声明不能吸收两个内部终止**；一条未绑定该位点的记录也不能验证别处的同类 stop；
     未验证、未解释的内部终止仍是 ERROR。
     registry 的最小 transl_except 记录：
     `{"gene":"cox1","codon":"TAA","amino_acid":"Trp","pos":"10..12",`
     `"transl_table":5,"taxon":"Lepidoptera","source":"DOI ...","rationale":"..."}`
     （或用 `"codon_index":4`，或 `"scope":"gene_wide"` 且不给位点）。

partial 的实现约束：

   - 代码直接读 Biopython 的 `BeforePosition`/`AfterPosition` **位置对象**（不解析 location 字符串），
     按链方向解释生物学的 5'/3' 端；**负链的 5' 端在高坐标**，跨原点 `join()` 的各段按转录顺序排列；
     位置对象**始终优先于**字符串兜底（仅当 parts 完全不带 fuzzy 信息时才用字符串），
     两者矛盾时以对象为准**并报** `PARTIAL_SOURCE_CONFLICT`，不静默取值；

## 9. COX1 BLAST 结果解析契约

查询片段为 400–5000 bp，结果只提供分类线索，不能单独确定物种。它请求 `FORMAT_TYPE=XML2`，并**同时**能读
XML2（`-outfmt 16`）与旧版 XML（`-outfmt 5`）两种方言；`FORMAT_OBJECT=SearchInfo` 按 NCBI 实际返回的
QBlast 文本（`RID =` / `Status=`）解析，XML 形式也兼容。**每个 `<Hsp>` 的 query/hit 坐标、identity、align-len、
bit-score 都必须存在且数值合法**：缺任一必需字段时整份结果按**格式故障**（退出码 3）处理 ——
缺少 subject 坐标就无法做方向/共线性/目标重用检查，不能靠“只剩单个 HSP”绕过验证。它报告每个候选的 **query coverage**（多 HSP 并集）与 identity，并把 RID 轮询与结果下载分开。
退出码把任务状态与判读状态分开：`0` 得判读 / `1` insufficient 或 no_match（分析结论）/ `2` 拒绝执行 /
`3` 网络或结果格式故障。多 HSP 的指标只有同时满足**不重叠**与**目标共线**才回总，检查分开做：
① 全部 HSP 必须同一目标链（真正的正负链混合属矛盾证据）；
② 按 `query_from` 排序后，`hit_from` 必须**沿目标链方向单调推进**（正链递增、**负链递减**）；
③ **query 侧与 target 侧都不能复用同一批碱基**：query 区间重叠记 `overlapping_hsps`，target 区间重叠记
`subject_overlap`（两条 query 分别打到目标同一段的重复/塌缩情形，仍不是两条独立比对）；
其中 **`max_span_ratio=3.0`（目标跨度 ≤ 3 倍比对长度）是额外的工程预警，不是 COX1 生物学标准，也不能代替方向检查**。
坐标约定经本地真实 `blastn` 实测锁定：连续负链 `q1..400→s1800..1401`、`q501..900→s1400..1001`（应**接受**）；
同两段但顺序倒置 `q1..400→s1400..1001`、`q500..900→s1801..1401`（应**拒绝**）。
回总失败时记 `overlapping_hsps` / `subject_overlap` / `non_collinear_hsps` + `CONFLICTING_ALIGNMENT`；若跳变同时触及目标两端，
另记 `CROSS_ORIGIN_CANDIDATE` —— 环状参考下可能是跨原点排列，但**必须**有明确坐标与结构证据，不得直接按连续线性比对接受。
这些情形都**不参与自动择优**，但 `AMBIGUOUS_ALIGNMENT` 的含义是“**多 HSP 指标无法可靠汇总**”，
**不是**“该 hit 不是有效候选”：CLI 会逐条打印被阻断候选的 HSP 明细（`q.. → s..  identity bits`），
并用 `--output-json` 把每个候选的全部 HSP 与 blocker 持久化；该文件无论判读结果如何都会写出。最优 accession 不等于已完成物种鉴定。
