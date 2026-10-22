# -*- coding: utf-8 -*-
"""libdock.fireprot —— FireProt 2.0 Web 服务客户端(蛋白热稳定性设计, 学术免费)
===================================================================================
FireProt (https://loschmidt.chemi.muni.cz/fireprotweb/) 是 Loschmidt 实验室的
全自动蛋白热稳定性设计平台。**无需账号/API key**, 只需提交 email + PDB ID(或自
定义结构), 服务端跑完整流程(同源建模 → 进化/能量/力场分析 → 突变建议)。

本模块用纯标准库 urllib 实现完整 API 流程(与 fetch.py / fold_api.py 同口径,
零第三方依赖):

  1. GET  /api/pdb/getJobId           -> 服务端分配 jobId
  2. GET  /api/pdb/processPdbId       -> 解析 PDB ID 得到生物单元/链
  3. GET  /api/CataliticResidues/getFastaSequence -> 提取 FASTA 序列
  4. POST /api/addNewJob              -> 提交 job(返回 jobId)
  5. GET  /api/jobInfo                -> 轮询状态
  6. GET  /api/getAllData             -> 结果 JSON(突变建议/稳定性预测)
  7. GET  /api/getZip                 -> 原始数据 zip(结果 PDB/报告)

状态机(job 整体): 1=New, 3=In progress, 5=Done, -1=Failed。
代理: 显式 proxy 参数 > DOCK_PROXY 环境变量 > 直连。

诚实边界: FireProt 是 Web 黑盒服务, 其打分口径(FoldX/Rosetta/进化共识等)与
DockBuddy 其它打分模型不可混比; 输出是「设计建议」而非绝对 ΔΔG。
"""
from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request
from typing import Any, Dict, List, Optional, Tuple

BASE_URL = "https://loschmidt.chemi.muni.cz/fireprotweb"

# FireProt 前端默认 advancedOptions(阈值; 与官方 UI 默认值一致)
DEFAULT_ADVANCED = {
    "maxNumSeqField": 200,
    "blastEvalField": 1e-10,
    "minIdentityField": 30,
    "maxIdentityField": 90,
    "majorityConsensusMinIdentity": 50,
    "ratioConsensusMinFreq": 40,
    "ratioConsensusMinRatio": 5,
    "globalConsensusGapThreshold": 50,
    "evolutionFoldXThreshold": 0.5,
    "energyFoldXThreshold": -1,
    "energyRosettaThreshold": -2,
    "multipointRosettaThreshold": -2,
    "limitForceFields": True,
    "zScoreThreshold": 3.5,
}

# FireProt 前端默认 selection(AUTOMATED 单点+多点; highRiskDesign 默认开)
DEFAULT_SELECTION = {
    "lowRiskDesign": False,
    "highRiskDesign": True,
    "ancestralDesign": False,
    "forceFieldMutations": True,
    "bFactorAnalysis": True,
    "consensusMutations": True,
    "isUserDefined": False,
}

# job 整体状态 -> 文本(与前端 Mr 枚举一致)
_STATUS_TEXT = {1: "New", 3: "In progress", 5: "Done", -1: "Failed"}


class FireProtError(RuntimeError):
    pass


def _default_proxy() -> str:
    return os.environ.get("DOCK_PROXY", "")


def _opener(proxy: str):
    handlers = []
    if proxy:
        handlers.append(urllib.request.ProxyHandler({"http": proxy, "https": proxy}))
    return urllib.request.build_opener(*handlers)


def _get(url: str, proxy: str, timeout: int) -> bytes:
    """GET 返回响应体字节。HTTP 错误抛 FireProtError。"""
    try:
        with _opener(proxy).open(url, timeout=timeout) as resp:
            return resp.read()
    except urllib.error.HTTPError as e:
        raise FireProtError(f"FireProt HTTP {e.code} ({url})") from e
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise FireProtError(f"FireProt 请求失败 ({url}): {e}") from e


def _post_json(url: str, payload: Dict[str, Any], proxy: str, timeout: int) -> bytes:
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=body, method="POST")
    req.add_header("Content-Type", "application/json")
    try:
        with _opener(proxy).open(req, timeout=timeout) as resp:
            return resp.read()
    except urllib.error.HTTPError as e:
        raise FireProtError(f"FireProt HTTP {e.code} ({url})") from e
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise FireProtError(f"FireProt 请求失败 ({url}): {e}") from e


def _get_text(url: str, proxy: str, timeout: int) -> str:
    return _get(url, proxy, timeout).decode("utf-8", "ignore").strip()


def _get_json(url: str, proxy: str, timeout: int) -> Dict[str, Any]:
    raw = _get(url, proxy, timeout).decode("utf-8", "ignore")
    try:
        return json.loads(raw)
    except json.JSONDecodeError as e:
        raise FireProtError(f"FireProt 返回非 JSON ({url}): {raw[:120]!r}") from e


def submit(pdb_id: str, email: str, *, job_title: str = "",
           selection: Optional[Dict[str, Any]] = None,
           advanced_options: Optional[Dict[str, Any]] = None,
           molecule: Optional[int] = None,
           proxy: str = "", timeout: int = 120) -> Dict[str, Any]:
    """提交 FireProt job(public PDB ID 输入), 返回 jobId 与输入摘要。

    参数:
      pdb_id: PDB ID(如 5XJH)。
      email: 收结果通知的邮箱(必填; FireProt 无账号, 仅用 email 限流/通知)。
      job_title: job 标题(缺省取 pdb_id)。
      selection: 分析类型开关(缺省 DEFAULT_SELECTION; AUTOMATED 单点+多点)。
      advanced_options: 阈值(缺省 DEFAULT_ADVANCED)。
      molecule: 生物单元编号(缺省取第一个 structures[0].molecule)。
      proxy: HTTP/HTTPS 代理(缺省 DOCK_PROXY 环境变量)。
      timeout: 单次请求超时(秒)。
    返回: {"job_id", "sequence", "chain_mapping", "pdb_id", "email"}。
    """
    p = proxy or _default_proxy()
    if not email or "@" not in email:
        raise ValueError("FireProt 需要有效 email(服务端强制校验)")
    title = (job_title or pdb_id).strip() or pdb_id

    # 1. 分配 jobId
    job_id = _get_text(f"{BASE_URL}/api/pdb/getJobId?jobId=", p, timeout)
    if not job_id:
        raise FireProtError("FireProt 未返回 jobId")

    # 2. 解析 PDB ID -> 生物单元/链
    st_json = _get_json(f"{BASE_URL}/api/pdb/processPdbId?pdbId={pdb_id}&jobId={job_id}",
                        p, timeout)
    structures = st_json.get("structures") or []
    if not structures:
        raise FireProtError(f"FireProt 未解析出 {pdb_id} 的结构(无效 PDB ID?)")
    st = structures[0]
    mol = molecule if molecule is not None else int(st.get("molecule", 1))
    chains = st.get("chains", "")
    original = st.get("original", "")
    chain_mapping = f"{chains}:{original}"

    # 3. 提取 FASTA 序列(指定链/生物单元)
    seq = _get_text(f"{BASE_URL}/api/CataliticResidues/getFastaSequence"
                    f"?jobId={job_id}&specifyChainsManually=false&selectedMM={mol}",
                    p, timeout)

    # 4. 提交 job
    body: Dict[str, Any] = {
        "jobId": job_id,
        "jobTitle": title,
        "email": email,
        "sequence": seq,
        "advancedOptions": dict(DEFAULT_ADVANCED if advanced_options is None
                                else advanced_options),
        "selection": dict(DEFAULT_SELECTION if selection is None else selection),
        "pdb": pdb_id,
        "model": False,          # Structure 输入(=0)
        "userDefined": [],
        "annotations": None,
        "uniprotId": "",
        "chainMapping": chain_mapping,
        "inputFileJson": json.dumps({"jobTitle": title}),
    }
    resp = _post_json(f"{BASE_URL}/api/addNewJob", body, p, timeout)
    out_job_id = resp.decode("utf-8", "ignore").strip()
    if not out_job_id:
        raise FireProtError("FireProt 提交被拒(可能超出并发 job 上限, 稍后再试)")
    return {"job_id": out_job_id, "sequence": seq, "chain_mapping": chain_mapping,
            "pdb_id": pdb_id, "email": email, "job_title": title}


def poll(job_id: str, proxy: str = "", timeout: int = 60) -> Dict[str, Any]:
    """查询 job 状态, 返回 jobInfo 原始 JSON(含 status / jobCalculations)。"""
    p = proxy or _default_proxy()
    info = _get_json(f"{BASE_URL}/api/jobInfo?jobId={job_id}", p, timeout)
    info["status_text"] = _STATUS_TEXT.get(info.get("status"), str(info.get("status")))
    return info


def wait(job_id: str, *, proxy: str = "", poll_interval: int = 30,
         max_wait: int = 7200, timeout: int = 60) -> Dict[str, Any]:
    """轮询直到 job 完成(status=5)或失败(status=-1)或超时。

    返回最终 jobInfo。注意 FireProt 完整流程可能需 30 分钟~数小时。
    """
    p = proxy or _default_proxy()
    deadline = time.time() + max_wait
    while time.time() < deadline:
        info = poll(job_id, p, timeout)
        status = info.get("status")
        if status == 5:
            return info
        if status == -1:
            raise FireProtError(f"FireProt job {job_id} 失败: {info}")
        time.sleep(poll_interval)
    raise FireProtError(f"FireProt job {job_id} 超时({max_wait}s 未完成), 可稍后用同一 jobId 继续 poll")


def get_results(job_id: str, proxy: str = "", timeout: int = 120) -> Dict[str, Any]:
    """获取完成的 job 结果 JSON(突变建议/稳定性预测等)。"""
    p = proxy or _default_proxy()
    return _get_json(f"{BASE_URL}/api/getAllData?jobId={job_id}", p, timeout)


def download_zip(job_id: str, out_path: str, proxy: str = "", timeout: int = 300) -> str:
    """下载 job 原始数据 zip(结果 PDB/报告), 写到 out_path。返回 out_path。"""
    p = proxy or _default_proxy()
    data = _get(f"{BASE_URL}/api/getZip?jobId={job_id}", p, timeout)
    if not data or len(data) < 22:  # zip 最小头 22 字节
        raise FireProtError(f"FireProt getZip 返回空/无效内容({len(data)} 字节)")
    with open(out_path, "wb") as f:
        f.write(data)
    return out_path


def submit_and_wait(pdb_id: str, email: str, *, job_title: str = "",
                    selection: Optional[Dict[str, Any]] = None,
                    advanced_options: Optional[Dict[str, Any]] = None,
                    out_dir: str = "", proxy: str = "",
                    poll_interval: int = 30, max_wait: int = 7200,
                    timeout: int = 120) -> Dict[str, Any]:
    """提交 + 轮询 + 拉取结果 + 下载 zip 的一站式流程。

    返回 {"job_id", "info", "results", "zip_path", ...}。结果 zip 落盘到 out_dir
    (缺省不落盘)。适合长任务后台跑。
    """
    sub = submit(pdb_id, email, job_title=job_title, selection=selection,
                 advanced_options=advanced_options, proxy=proxy, timeout=timeout)
    job_id = sub["job_id"]
    info = wait(job_id, proxy=proxy, poll_interval=poll_interval,
                max_wait=max_wait, timeout=timeout)
    results = get_results(job_id, proxy=proxy, timeout=timeout)
    out: Dict[str, Any] = {"job_id": job_id, "info": info, "results": results,
                           "zip_path": "", **sub}
    if out_dir:
        import os as _os
        _os.makedirs(out_dir, exist_ok=True)
        out["zip_path"] = download_zip(job_id, _os.path.join(out_dir, f"{job_id}.zip"),
                                       proxy=proxy, timeout=timeout)
    return out


__all__ = ["FireProtError", "BASE_URL", "DEFAULT_ADVANCED", "DEFAULT_SELECTION",
           "submit", "poll", "wait", "get_results", "download_zip", "submit_and_wait"]
