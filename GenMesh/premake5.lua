newoption {
   trigger = "openvdb-path",
   description = "Path to OpenVDB 4.0+ library"
}

newoption {
   trigger = "tbb-path",
   description = "Path to Intel TBB library"
}

sources = {
   "src/main.cpp",
}

-- premake5.lua
workspace "GenMeshWorkspace"
   configurations { "Release", "Debug" }

   platforms { "native", "x64", "x32" }

   -- Use c++11
   flags { "c++11" }

   includedirs { "./include" }

   includedirs { "./src" }

   -- A project defines one build target
   project "GenMesh"
      kind "ConsoleApp"
      language "C++"
      files { sources }

      if os.is("Linux") then
         -- TBB
         if _OPTIONS["tbb-path"] then
            includedirs { _OPTIONS["tbb-path"] .. "include" }
            libdirs { _OPTIONS["tbb-path"] .. "lib" }
         end
        
         -- OpenVDB
         includedirs { _OPTIONS["openvdb-path"] .. "/include/" }
         libdirs { _OPTIONS["openvdb-path"] .. "/lib/" }
         links { "openvdb" }

         links {"pthread", "dl", "Half", "tbb", "boost_iostreams"}
      end

      configuration "Debug"
         defines { "DEBUG" } -- -DDEBUG
         symbols "On"
         targetname "genMesh_debug"

      configuration "Release"
         -- defines { "NDEBUG" } -- -NDEBUG
         symbols "On"
         flags { "Optimize" }
         targetname "genMesh"
