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
   而**格式非法（§2.1–2.8）导致加载失败并退出 1**，不进入 REVIEW，也不做部分生效。

## 3. 退出码语义（不得误读）

| 工具 | 码 | 含义 | **不是** |
|---|---|---|---|
| `annot_check.py` | `0` / `1` / `2` | 无发现 / 有错误（含参数错误）/ 仅待核查 | `2` **不是**通过；`case-validate` 的通过也不是科学验收 |
| `cox1_id.py` | `0` / `1` / `2` / `3` | 得判读 / `insufficient` 或 `no_match` / 拒绝执行 / 网络或结果格式故障 | `3` **不是**"没找到相似序列" |
| `blast_genes.py` | `0` / `1` | 全部基因唯一命中 / 任一基因无唯一命中或命中歧义（**不写部分结果**） | `1` 不是"基因真的缺失" |
| `mitos2_to_genbank.py` | `0` / 非 0 | 写出 `<out.gb>` / 缺文件、坐标越界、`result.fas` 与坐标不一致、未知 feature 类型 | 非 0 时**不得**留下半成品 |
| `depth_analysis.py` | `0` / `1` / `2` | 正常 / 输入或 BAM 错误 / 存在低覆盖区 | `2` 不是"组装错误" |
| `circularize.py` | `2` / `1` | REVIEW/需人工 / 接缝证据不足 | 二者都**不是**"已确认环化" |

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
