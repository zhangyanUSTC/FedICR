# 重建 PISA 2022 数据

[English](README.md)

本项目提供从 OECD 官方学生问卷文件重建分析数据的代码，不随仓库分发原始或处理后的学生级数据。模拟实验不需要这些文件。

## 1. 下载并放置源数据

打开 [OECD PISA 数据下载页](https://webfs.oecd.org/pisa2022/index.html)，在 **PISA 2022 Data → SAS Data Files (Compressed)** 下下载 **Student questionnaire data file**。解压后，将 SAS 文件放在：

~~~text
data/CY08MSP_STU_QQQ.SAS7BDAT
~~~

请选择 2022 年的学生问卷文件，不要选择其他年份、学校、教师或认知测验题目文件。源文件有数 GB，需要预留足够的内存和磁盘空间。[OECD PISA 2022 数据库页面](https://www.oecd.org/en/data/datasets/pisa-2022-database.html)也提供下载入口。

## 2. 安装依赖并重建数据

在包含 run.R 的项目根目录运行：

~~~sh
Rscript scripts/install_packages.R
Rscript scripts/pisa_data_rebuild.R
Rscript pisa/prepare_pv1.R
~~~

第一条数据准备命令生成十套 PV 数据，第二条从中提取 PV1，不改变样本。

| 生成文件 | 内容 | 使用任务 |
| --- | --- | --- |
| data/pisa_federated_data_pv.rds | 十套匹配的 PV，每个客户端 35 列 | pisa-sensitivity |
| data/pisa_federated_data_pv1.rds | 仅 PV1，每个客户端 8 列 | pisa-cva、pisa-compare |

两个准备脚本都不会覆盖已存在的输出。如果已有相应的处理后数据文件，可跳过该准备步骤。

## 3. 重建过程

脚本读取 CNT、ANXMAT、MISCED、FISCED 以及阅读、科学、数学的全部 30 个 PV 列；统一删除这些字段中有缺失值的记录，添加名为 1 的截距列，再按 CNT 拆分为国家或经济体的 data.frame 列表。客户端内部保留学生的原始顺序。不对 PV 求平均，也不标准化所保留的变量；分析不使用最终学生权重或重复权重。

CNT 标识客户端。第 j 套 PV 分析以数学为响应变量，以阅读、科学、ANXMAT、MISCED、FISCED 和截距为预测变量。请保留生成的列顺序及各套 PV 的共同学生样本，不要分别筛选或排序不同的 PV 数据。

## 4. 检查样本并运行 PISA

重建完成后会打印验证摘要。论文使用的样本为 **75 个客户端、451,851 名学生、35 列、10 套 PV**。PV1 保留相同的客户端和学生，共 8 列。

如果数量不同，应核实源文件版本和缺失值处理，不要为了达到预期数量而增删记录。

~~~sh
Rscript run.R --task pisa-all --workers 10
~~~

若处理后数据存放在其他目录，请添加 --data 并指定该目录，包含空格的路径需加引号。各任务的运行方式及结果解读见[主复现说明](../README.md)。
