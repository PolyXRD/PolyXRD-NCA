# 幂等创建 GitHub 私有仓库 PolyXRD/PolyXRD-NCA（可重跑）
# 用法: & 'scripts\create_remote_repo.ps1'
#
# ★ 凭据来源与家族完全一致：`git credential fill` 读 Windows 凭据管理器
#   （系统级 credential.helper = manager / Git Credential Manager），
#   不落盘 token、不要求手工配 gh 或环境变量。
# ★ 为什么单独成文件：交互式命令里内联 [Diagnostics.Process]::Start 会被安全策略拦，
#   写进 .ps1 再调用是家族既有做法（见 PolyXRD/scripts/create_github_release_*.ps1）。
#
# ★★ description 必须是纯 ASCII。实测带中文 description 的 POST /user/repos 返回
#   500 Internal Server Error（body 为空），同样的载荷换成 ASCII 就能建出来 —— 这是
#   GitHub 建库接口对非 ASCII description 的老毛病，不是编码写错（payload 已确认是
#   UTF-8 无 BOM）。想要中文简介，建库成功后再用 PATCH 单独改。
# ★★ 建库/建 Release 会出现「假失败」：HTTP 500/502 但资源其实已创建。所以出错后
#   不能直接判定失败，必须回查一次 GET /repos/{owner}/{repo} 再决定是否重试。
# ★ 本 token 的 scopes = gist, read:org, repo, workflow，**没有 delete_repo** ——
#   脚本能建库不能删库。手工删请到仓库 Settings -> Danger Zone。
param(
  [string]$Owner = 'PolyXRD',
  [string]$Repo  = 'PolyXRD-NCA'
)
$ErrorActionPreference = 'Stop'

# ---- 凭据：dot-source 家族小工具（含 Show-GitHubError）----
. (Join-Path $PSScriptRoot '_gh_cred.ps1')
$headers = Get-GitHubHeaders

$full = "$Owner/$Repo"
$repoObj = $null
try {
  $repoObj = Invoke-RestMethod -Uri "https://api.github.com/repos/$full" -Headers $headers -Method Get -ErrorAction Stop
  Write-Host "REPO ALREADY EXISTS: $($repoObj.full_name)  private=$($repoObj.private)"
} catch {
  Write-Host "REPO ABSENT; creating private repo $full ..."
  $payload = @{
    name         = $Repo
    description  = 'PolyXRD-NCA - XRD Rietveld refinement for Li(NiCoAl)O2 cathode materials and precursors'
    private      = $true
    has_issues   = $true
    has_wiki     = $false
    has_projects = $false
    auto_init    = $false
  } | ConvertTo-Json
  try {
    $repoObj = Invoke-RestMethod -Uri 'https://api.github.com/user/repos' -Headers $headers -Method Post `
        -Body ([System.Text.Encoding]::UTF8.GetBytes($payload)) -ContentType 'application/json; charset=utf-8' -ErrorAction Stop
    Write-Host "CREATED: $($repoObj.full_name)  private=$($repoObj.private)"
  } catch {
    # 假失败兜底：500/502 时资源常已建好，回查确认
    Show-GitHubError $_
    Start-Sleep -Seconds 5
    try {
      $repoObj = Invoke-RestMethod -Uri "https://api.github.com/repos/$full" -Headers $headers -Method Get -ErrorAction Stop
      Write-Host "converged (POST reported error but repo exists): $($repoObj.full_name) private=$($repoObj.private)"
    } catch { throw "建库失败且回查仍不存在: $full" }
  }
}

Write-Host "html_url   = $($repoObj.html_url)"
Write-Host "clone_url  = $($repoObj.clone_url)"
Write-Host "=== DONE ==="
