# 任务证据记录（版本 3）

用于本次诊断的输入追溯、执行记录、逐条结论和报告。由代理维护 JSON，用户阅读报告即可。
程序不执行登记的命令，也不认证命令、软件版本或复核人身份；哈希只标识文件内容。
不包含自动修复、跨任务学习或旧案例迁移。适用于 Linux/POSIX 本地任务目录。

## 操作顺序

以下命令在 skill 根目录运行，路径换成实际路径。所有分析产物放在用户任务目录，
原始文件只读；保留已登记文件，后续分析使用不同输出文件名。

```bash
python3 scripts/case.py init /work/sample-review --case-id sample-review \
  --issue 'cox1 边界待核查' --taxon Insecta \
  --input assembly /data/assembly.fa --input annotation /data/annotation.gb
python3 scripts/case.py event /work/sample-review --record /work/event.json
python3 scripts/case.py conclusion /work/sample-review --record /work/conclusion.json
python3 scripts/case.py validate /work/sample-review --verify-files
python3 scripts/case.py report /work/sample-review
```

依赖 `jsonschema>=4.18`，安装到运行脚本的 Python 环境。退出码 0 表示记录操作成功，
1 表示参数、格式、关联或文件校验失败，3 表示缺依赖。成功不等于科学验收。
`init` 输出输入 ID（I1、I2…）和哈希；事件引用实际输出的 ID。
新增输入：`python3 scripts/case.py input /work/sample-review --role reads --path /data/reads.fq`。
role 可描述用途；reads 支持检查识别 `reads`、`bam`，竞争参考使用 `competitive_reference`。
遗传密码表及其依据、参考数据库版本等上下文写入相应事件 summary。

先按当前异常实际运行必要检查，再提交事件；不要将下面的示例当作已经执行的事实。
事件中的退出码是分析工具的退出码，完成不要求为 0，例如 annot_check 的 2 表示待核查。
`failed` 表示工具执行失败，`not_run` 的退出码必须为 null；两者只能支持 tool_execution 类结论。

## 事件 JSON 示例

将工作目录、版本、日志路径和观察替换为真实值；command 是参数数组，不是 shell 字符串。
重定向日志的操作应单独如实记录，不能把 `>` 当作分析脚本的参数。

```json
{
  "id": "E1",
  "action": "annotation_check",
  "command": ["python3", "scripts/annot_check.py", "/data/annotation.gb", "--table", "5"],
  "working_directory": "/path/to/mito-assembly-triage",
  "tool_version": "填写实际脚本版本或提交号",
  "exit_code": 2,
  "execution_status": "completed",
  "input_ids": ["I2"],
  "summary": "填写实际观察；说明密码表确认依据及竞争解释",
  "artifacts": [{"id": "check_log", "role": "log", "path": "/work/sample-review/check.log"}]
}
```

产物仅提交 id/role/path，程序计算大小和 SHA-256，登记时自动写 recorded_at。
artifact ID 在所属事件内唯一，event ID 在整个任务内唯一；已有事件不能覆盖。
completed 事件必须有已登记输入及真实产物。证据引用精确到 event_id + artifact_id。

## 结论 JSON 示例

下面展示有部分支持但仍未解决的写法；observation 必须与实际产物相符。

```json
{
  "id": "C1",
  "claim": "cox1 边界仍需进一步确认",
  "conclusion_type": "annotation_quality",
  "status": "UNRESOLVED",
  "confidence": "moderate",
  "reads_support": "NOT_ASSESSED",
  "evidence_for": [{"event_id": "E1", "artifact_id": "check_log", "observation": "填写支持这一命题的实际观察"}],
  "evidence_against": [],
  "not_tested": ["尚未比较近缘参考的保守蛋白边界"],
  "scope": "当前 cox1 注释边界",
  "limitations": ["此检查不验证样本碱基真实性"],
  "rationale": "说明为何现有证据尚不能区分竞争边界",
  "next_steps": ["比较可靠近缘参考的蛋白边界"],
  "decided_by": "填写实际判定人或代理标识",
  "review": {"level": "ASSIST", "status": "pending", "reviewer": null, "note": "等待边界复核"},
  "competitive_reference_ids": []
}
```

conclusion_type 支持 assembly_quality、annotation_quality、gene_identity、sequence_accuracy、
biological_interpretation、tool_execution。confidence 支持 high、moderate、low、not_assessable。
每个命题单独判定；RESOLVED 和 high 必须有支持证据；RESOLVED 不得配 not_assessable。
NO_CHANGE 只表示不修改，缺证据时说明理由；UNRESOLVED 可包含部分支持及冲突。
review.level 支持 AUTO/ASSIST/EXPERT；status 支持 not_requested/pending/reviewed/disputed，
后两者需要 reviewer。复核状态不替代科学状态。

READS_CONSISTENT 必须引用使用 reads/bam 输入的支持事件。READS_DISCRIMINATING 还要求
登记竞争参考 ID 并关联用于比较的支持事件。角色及关联检查只是必要条件，不能证明
实际完成竞争比对或排除 NUMT；须核对真实命令、比对结果和适用范围。
正常案例 normal_validation_case 不允许 UNRESOLVED；诊断期间默认 abnormal_case，
不要为得到“正常”标签提前做决定。工具故障可使用 tool_failure_case。

## 修订与文件完整性

更新已有结论用 `conclusion ... --record ... --update`，旧结论进入 conclusion_history。
input/event/conclusion 均可加 `--expected-revision N` 防止提交过期决定。
case.json 原子保存所有事件与结论，进程锁避免合作进程丢失更新；不要手工改这个文件。
它不是签名审计系统，不能防止有写权限者直接篡改历史。

所有写入和 report 都重新核对已登记文件；文件缺失或变化会阻止操作。
不要覆盖旧产物；保存原始版本，再将新版本作为新输入/事件登记。
文件很大时哈希校验需要相应 I/O 时间。validate 默认仅校验格式及关联，
`--verify-files` 才核对文件；移动任务时须保留其绝对路径或另建任务重新登记。

默认报告是任务目录 report.md，可用 `--output` 指定位置；禁止覆盖输入、产物或账本。
报告标明修订号、账本哈希、每条证据及未测项。失败时原有报告可能仍存在，
必须检查退出码与修订号，不得将旧报告当作本次成功输出。
历史 schemas/case.schema.json 与 experience.py 文档不是本接口；旧格式拒绝隐式迁移。

## 修复候选与修改后验证

具体候选文件由现有诊断/转换工具生成，必须写入新路径。完成生成后，先将实际命令及
候选文件作为产物登记到 `event`，再登记候选元数据：

```json
{
  "id": "R1",
  "change_type": "annotation",
  "base_input_id": "I2",
  "producer_event_id": "E2",
  "artifact_id": "candidate_gb",
  "rationale": "实际证据支持调整 cox1 边界",
  "expected_effect": "使翻译边界与独立近缘参考一致",
  "risks": ["需要确认类群遗传密码表", "不能验证样本碱基"]
}
```

```bash
python3 scripts/case.py candidate /work/sample-review --record /work/candidate.json
```

账本会校验生成事件确实引用基础输入，候选文件与原输入分离且内容不同，并将候选登记为
新输入 ID（例如 I3）。随后按候选类型执行**独立的**必要检查；检查事件的 `input_ids`
必须含 I3，并如实登记输出、退出码和工具版本。可将精确的统一文本差异也作为该事件产物
保存，以便审阅；差异文本本身不代替生物学验证。

验证后提交记录：

```json
{
  "candidate_id": "R1",
  "event_id": "E3",
  "result": "pass",
  "criteria": "预先定义：显式密码表下无非预期内部终止，边界与两个独立近缘参考一致",
  "notes": "写明实际检查结果及限制"
}
```

```bash
python3 scripts/case.py candidate-check /work/sample-review --record /work/check.json
python3 scripts/case.py report /work/sample-review
```

`result` 为 `pass`、`fail` 或 `inconclusive`。账本按记录结果计算候选状态：无验证为
`proposed`；所有检查均 pass 才为 `verified`；所有检查均 fail 才为 `rejected`；
其他混合或未定结果均为 `inconclusive`。这里的 `verified` 只表示登记事件报告通过所述标准，不能认证事件真实性、
标准是否充分或生物学结论正确。报告会显示标准和事件，需由用户/专家审阅。

注释候选须使用明确确认的密码表并跑独立注释/翻译检查；边界参考证据不验证样本碱基。
碱基替换须有样本 reads 支持，并检查碱基质量、MAPQ、独立模板、链向、竞争比对与混合信号。
结构候选须根据读长与重复特性验证每个新增接缝和闭合连接；候选组装图不证明物理结构。
`sequence` 或 `structure` 候选在缺少相应原始证据时可以保留为 proposed/inconclusive，不能升级。
本流程记录候选与验证，不会采纳候选、覆盖原件或发布结果。

## 按清单生成注释或碱基候选

`scripts/apply_candidate.py` 只执行用户明确列出的改动，不从同源参考自动推断碱基、边界或注释字段。
每项都要写证据引用（建议 `E1/artifact_id`）和理由；源值不匹配、坐标重复或输出覆盖
输入或已存在的候选/manifest 时拒绝写入。工具额外输出包含源/编辑清单/候选 SHA-256 的 JSON manifest，须与候选
一起登记为生成事件产物。manifest 说明 `candidate_only_not_validated`。

碱基 edits 示例（位置以相应 FASTA record 为准，1-based）：

```json
{"changes":[{"sequence_id":"contig_mt","position":417,"old":"G","new":"A","evidence":"E2/read_support","rationale":"覆盖、质量和独立模板均支持 A；竞争比对已检查"}]}
```

```bash
python3 scripts/apply_candidate.py base assembly.fa edits.json candidate.fa --manifest candidate-changes.json
```

只接受新碱基为明确 A/C/G/T 的替换，旧值可为 IUPAC 模糊碱基；多记录 FASTA 必须给出 `sequence_id`。工具会核对每个旧碱基，
保留其他记录与模糊字符。它不会分析 reads；变更依据必须来自任务中已经检查过的证据。
碱基候选要登记为 `change_type: sequence`，只有引用 reads 或 BAM 的验证事件才能记 pass。
验证时按该命题核查碱基质量、MAPQ、独立模板、链向、混合信号和可用的竞争参考。
缺核参考时说明 NUMT 排除限制；不能把覆盖或单参考比对当作充分证据。

CDS 边界 edits 示例，GenBank feature 通过唯一 locus_tag 选择：

```json
{"changes":[{"selector":{"qualifier":"locus_tag","value":"cox1"},"old":{"start":1,"end":1542,"strand":1},"new":{"start":1,"end":1539,"strand":1},"old_codon_start":1,"new_codon_start":1,"evidence":"E3/homology","rationale":"两个独立近缘蛋白比对共同支持终止边界"}]}
```

```bash
python3 scripts/apply_candidate.py cds-boundary annotation.gb edits.json candidate.gb --table 5 --manifest candidate-changes.json
```

坐标为 1-based inclusive；密码表必须按类群确认并显式提供。旧、新 codon_start 都必须显式列出（值 1/2/3），旧值须匹配 GenBank；改起点时需依据证据选择新值。工具只改唯一命中的简单、
精确 CDS 的位置，不编辑 fuzzy/compound location；记录翻译前后序列供审查，不据此自动
接受或判断非典型例外。随后重新运行 `annot_check.py <candidate.gb> --table <确认表>`，
将候选、变更 manifest、检查日志分别登记；需保留 partial、例外和注释完整性的人工判断。

GenBank 来源、名称、标识符等 qualifier 可逐项按旧值精确编辑。`source` 选择器为
`{"type":"source"}`；其他特征用 `type` 加唯一的 `gene`、`locus_tag` 或 `product` 选取。
`old`/`new` 是字符串数组，空 `new` 删除字段；同一特征同一字段不能重复编辑。
修改 source `/organism` 时工具同步 GenBank 顶层 `SOURCE` 和 `ORGANISM`。例如：

```json
{"changes":[{"selector":{"type":"source"},"qualifier":"organism","old":["参考物种"],"new":["样本物种"],"evidence":"用户样本鉴定与原始来源记录","rationale":"移除参考物种身份"}]}
```

```bash
python3 scripts/apply_candidate.py qualifier annotation.gb edits.json candidate.gb --manifest candidate-changes.json
```

只改来源字段不能证明物种鉴定；从参考转来的采集地、日期、单倍型和蛋白 accession
也须逐项核对。生成候选后逐碱基比对 FASTA、确认全部基因坐标未变，并复跑注释质检。

结构候选使用 `circularize.py` 等现有结构工具产出后登记为 `change_type: structure`。
每个新增连接与最终尾首闭合都要单独验证。`candidate-check pass` 对碱基和结构候选要求
关联 reads/BAM 输入；序列/结构工具退出码非零不能记 pass。是否检查了所有必要位点、
竞争解释与验收标准是否充分，仍须由人审阅记录和原始产物。
