# 可克隆源码交付

目标 GitHub 仓库包含本次源代码、依赖锁、必要案例和历史版 ANDES 许可/源码：

```bash
git clone https://github.com/HardyHaoranLiu/CloudDRTS.git
cd CloudDRTS
git log -1
```

完整程序包 `CloudDRTS_Experiments_v7.zip` 中的 `CloudDRTS_Source_v7.bundle` 是可克隆的离线 v7 源码副本。大体积的生成结果在完整程序包中；Git 仓库通过入口重新生成。离线克隆方式：

```bash
git clone CloudDRTS_Source_v7.bundle CloudDRTS_Source_v7
cd CloudDRTS_Source_v7
git log -1
```

先按 README 建立自编环境并运行 `CloudDRTS/run_all.py`，再建立独立 ANDES 环境运行 `external_run/run_reproduction.py`。生成阶段说明时需要两边结果及 `requirements_notes.txt`；阶段一中英文说明另需 XeLaTeX，命令见 README。

所有代码以相对路径定位。离线 v7 源码提交为 `bc1fe4e5ce3dcb9f4f62cd60a8105631b24b648b`；完整程序包中的 `GIT_SOURCE.json` 记录提交号与逐文件校验值。Git bundle 完整性已经通过 `git bundle verify` 和克隆后逐文件 SHA-256 核对。

GitHub 版本保留该 v7 源码提交，并补充在线克隆说明。运行结果及原始实验状态以 README 中声明的适用范围和未完成项为准。
