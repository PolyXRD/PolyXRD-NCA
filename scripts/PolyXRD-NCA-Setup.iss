; ============================================================
; PolyXRD-NCA Setup - Inno Setup Script
; 使用: iscc /DAppVersion=0.1.0 /O"installer_output" /F"PolyXRD-NCA-v0.1.0-Setup" scripts\PolyXRD-NCA-Setup.iss
;       或直接跑 build_setup.bat（自动读 nca/version.py 的版本号并收尾 dist）。
;
; 语言: 简体中文 / 繁體中文 / English / 日本語
;   ⚠ Inno Setup 官方发行包【不自带】中文 .isl (仅 29 种语言),
;     故中文语言文件沿用家族仓库 scripts\languages\ 下的副本。
; ============================================================

#ifndef AppVersion
  #define AppVersion "0.1.0"
#endif

#define AppName "PolyXRD NCA"
#define AppPublisher "PolyXRD Team"
#define AppContact "sshztx@outlook.com"
#define AppExeName "PolyXRD-NCA.exe"
;# 家族的「一个 exe + 一个 Setup」口径（T8.5 对齐）：
;#   打包 payload 是 build_exe.bat 产出的**单文件**内部评估版，
;#   安装到 {app} 时改名为 PolyXRD-NCA.exe（用户看到的名字不带版本号）。
;#   历史上这里是 onedir 目录 ..\dist\PolyXRD-NCA，Setup 需递归拷 500+ 文件，
;#   改成单文件后 Setup 只装 1 个 exe，安装/卸载都干净。
#define BuildExe "..\dist\PolyXRD-NCA-Internal-v" + AppVersion + ".exe"
;# 中文语言文件所在目录 (相对本 .iss 文件) —— 自家族仓库
;#   D:\Project\XRD\PolyXRD\scripts\languages\ 复制而来的副本，
;#   使本仓库打包自包含 (不依赖兄弟目录存在)。
#define LangDir "languages"

[Setup]
AppId={{34FCCE89-4065-42C8-BF5D-AB40CCC09F1C}
; AppId 必须与兄弟项目不同 —— 从 PolyXRD-NCM 拷贝模板时曾原样带过来，
; 后果是装 NCA 会被 Inno 当成「升级 NCM」而覆盖/卸载 NCM。本项目的
; AppId 于 v0.1.0 重新生成，禁止再从别的项目复制本行。
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher={#AppPublisher}
AppPublisherURL=mailto:{#AppContact}
AppSupportURL=mailto:{#AppContact}
AppUpdatesURL=https://github.com/PolyXRD
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
AllowNoIcons=yes
OutputDir=..\installer_output
OutputBaseFilename=PolyXRD-NCA-v{#AppVersion}-Setup
SetupIconFile=..\nca\resources\app-icon.ico
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
; lowest = 默认按当前用户安装 (无需管理员); 允许用户在向导里切换为
;   "为所有用户安装", 避免装到需要管理员权限的位置时权限不足。
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
DisableProgramGroupPage=yes
DisableDirPage=no
UninstallDisplayIcon={app}\{#AppExeName}
VersionInfoVersion={#AppVersion}
VersionInfoDescription={#AppName} v{#AppVersion} Setup
VersionInfoCompany={#AppPublisher}
VersionInfoCopyright=(C) 2026 PolyXRD Team
SetupLogging=yes

[Languages]
; 中文放最前: 系统语言为中文时向导自动选中中文
Name: "chinesesimplified"; MessagesFile: "{#LangDir}\ChineseSimplified.isl"
Name: "chinesetraditional"; MessagesFile: "{#LangDir}\ChineseTraditional.isl"
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "japanese"; MessagesFile: "compiler:Languages\Japanese.isl"

[CustomMessages]
chinesesimplified.AdditionalIcons=附加图标:
chinesesimplified.CreateDesktopIcon=创建桌面快捷方式(&D)
chinesesimplified.CreateQuickLaunchIcon=创建任务栏快捷方式(&Q)
chinesesimplified.LaunchApp=启动 {#AppName}
chinesesimplified.UninstallShortcut=卸载 {#AppName}
chinesetraditional.AdditionalIcons=附加圖示:
chinesetraditional.CreateDesktopIcon=建立桌面捷徑(&D)
chinesetraditional.CreateQuickLaunchIcon=建立工作列捷徑(&Q)
chinesetraditional.LaunchApp=啟動 {#AppName}
chinesetraditional.UninstallShortcut=反安裝 {#AppName}
english.AdditionalIcons=Additional Icons:
english.CreateDesktopIcon=Create a &desktop shortcut
english.CreateQuickLaunchIcon=Create a &Quick Launch shortcut
english.LaunchApp=Launch {#AppName}
english.UninstallShortcut=Uninstall {#AppName}
japanese.AdditionalIcons=追加アイコン:
japanese.CreateDesktopIcon=デスクトップショートカットを作成(&D)
japanese.CreateQuickLaunchIcon=クイックランチショートカットを作成(&Q)
japanese.LaunchApp={#AppName} を起動
japanese.UninstallShortcut={#AppName} をアンインストール

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "{#BuildExe}"; DestDir: "{app}"; DestName: "{#AppExeName}"; Flags: ignoreversion

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\{#AppExeName}"; IconFilename: "{app}\{#AppExeName}"
Name: "{group}\{cm:UninstallShortcut}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExeName}"; IconFilename: "{app}\{#AppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExeName}"; Description: "{cm:LaunchApp}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
Type: filesandordirs; Name: "{localappdata}\PolyXRD-NCA"; Check: not DirExists('{userappdata}\PolyXRD-NCA')
Type: filesandordirs; Name: "{userappdata}\PolyXRD-NCA"
