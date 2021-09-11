workspace "IGameKernel" -- 解决方案
	architecture "x64"
    configurations
    {
        "Debug",
        "Release"
    }
	location "build"

project "Mesh" -- 项目
    kind "ConsoleApp" -- 控制台应用
    language "C++"
	targetdir "bin/Mesh/%{cfg.buildcfg}"
	location "build/Mesh"
	includedirs "./Mesh"
    files 
    {
		"./Mesh/Kernel/**.h",
        "./Mesh/Kernel/**.cpp", 
		"./Mesh/IO/**.h",
        "./Mesh/IO/**.cpp",
		"./Mesh/Tools/**.h",
        "./Mesh/Tools/**.cpp",
		"./Mesh/Application/**.h",
        "./Mesh/Application/**.cpp",
    }
filter "configurations:Debug"
      defines { "DEBUG" }
      symbols "On"

   filter "configurations:Release"
      defines { "NDEBUG" }
      optimize "On"