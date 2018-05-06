sources = {
   "src/main.cpp",
}

-- premake5.lua
workspace "GenTexWorkspace"
   configurations { "Release", "Debug" }

   platforms { "native", "x64", "x32" }

   -- Use c++11
   flags { "c++11" }

   includedirs { "./include" }

   includedirs { "./src" }

   -- A project defines one build target
   project "GenTex"
      kind "ConsoleApp"
      language "C++"
      files { sources }

      if os.is("Linux") then
         links { "GL", "GLEW", "glfw" }
      end

      configuration "Debug"
         defines { "DEBUG" } -- -DDEBUG
         symbols "On"
         targetname "genTex_debug"

      configuration "Release"
         -- defines { "NDEBUG" } -- -NDEBUG
         symbols "On"
         flags { "Optimize" }
         targetname "genTex"
