-- premake5.lua — minimal, system-installed OpenVDB + TBB only

workspace "GenMesh"
  configurations { "Debug", "Release" }
  location "build"

project "GenMesh"
  kind "ConsoleApp"
  language "C++"
  cppdialect "C++17"            -- 必要なら "C++20" に変更可
  targetdir ("%{wks.location}/%{cfg.buildcfg}/bin")
  objdir    ("%{wks.location}/%{cfg.buildcfg}/obj")

  files {
    "src/**.cpp",
    "include/**.hpp"
  }
  includedirs { "include" }

  filter "system:linux"
    links { "openvdb", "tbb" }  -- 余計な手動リンクはしない（最小）
    buildoptions { "-Wall", "-Wextra" }
    defines { "TINYGLTF_NO_STB_IMAGE", "TINYGLTF_NO_STB_IMAGE_WRITE" }

    -- 実行時に /usr/local/lib などが必要なら rpath を“任意で”付与:
    -- linkoptions { "-Wl,-rpath,/usr/local/lib" }

  filter "configurations:Debug"
    symbols "On"
    defines { "DEBUG" }

  filter "configurations:Release"
    optimize "Speed"
    defines { "NDEBUG" }

  filter {}
