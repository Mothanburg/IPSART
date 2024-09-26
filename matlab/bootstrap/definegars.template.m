%% About definegars.m
% This file defines the MATLAB interface to the library |gars|.
%
% Commented sections represent C++ functionality that MATLAB cannot automatically define. To include
% functionality, uncomment a section and provide values for ["height" "width"], <DIRECTION>, etc. For more
% information, see helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') to "Define MATLAB Interface for C++ Library".



%% Setup
% Do not edit this setup section.
function libDef = definegars()
libDef = clibgen.LibraryDefinition("garsData.xml");

%% OutputFolder and Libraries 
libDef.OutputFolder = "H:\Code\GaRS\native\build";
libDef.Libraries = "H:\Code\GaRS\native\build\windows\x64\release\GaRS.lib";

%% C++ function |GaRSTestOpenCL| with MATLAB name |clib.gars.GaRSTestOpenCL|
% C++ Signature: int GaRSTestOpenCL()

GaRSTestOpenCLDefinition = addFunction(libDef, ...
    "int GaRSTestOpenCL()", ...
    "MATLABName", "clib.gars.GaRSTestOpenCL", ...
    "Description", "clib.gars.GaRSTestOpenCL Representation of C++ function GaRSTestOpenCL."); % Modify help description values as needed.
defineOutput(GaRSTestOpenCLDefinition, "RetVal", "int32");
validate(GaRSTestOpenCLDefinition);

%% C++ function |RefinedLeeFilter3x3| with MATLAB name |clib.gars.RefinedLeeFilter3x3|
% C++ Signature: int RefinedLeeFilter3x3(long nLooks,long height,long width,float const * c11,float const * c22,float const * c33,float const * c12_r,float const * c13_r,float const * c23_r,float const * c12_i,float const * c13_i,float const * c23_i,float * outC11,float * outC22,float * outC33,float * outC12_r,float * outC13_r,float * outC23_r,float * outC12_i,float * outC13_i,float * outC23_i)

RefinedLeeFilter3x3Definition = addFunction(libDef, ...
   "int RefinedLeeFilter3x3(long nLooks,long height,long width,float const * c11,float const * c22,float const * c33,float const * c12_r,float const * c13_r,float const * c23_r,float const * c12_i,float const * c13_i,float const * c23_i,float * outC11,float * outC22,float * outC33,float * outC12_r,float * outC13_r,float * outC23_r,float * outC12_i,float * outC13_i,float * outC23_i)", ...
   "MATLABName", "clib.gars.RefinedLeeFilter3x3", ...
   "Description", "clib.gars.RefinedLeeFilter3x3 Representation of C++ function RefinedLeeFilter3x3."); % Modify help description values as needed.
defineArgument(RefinedLeeFilter3x3Definition, "nLooks", "int32");
defineArgument(RefinedLeeFilter3x3Definition, "height", "int32");
defineArgument(RefinedLeeFilter3x3Definition, "width", "int32");
defineArgument(RefinedLeeFilter3x3Definition, "c11", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "c22", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "c33", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "c12_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "c13_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "c23_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "c12_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "c13_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "c23_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC11", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC22", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC33", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC12_r", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC13_r", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC23_r", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC12_i", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC13_i", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter3x3Definition, "outC23_i", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineOutput(RefinedLeeFilter3x3Definition, "RetVal", "int32");
validate(RefinedLeeFilter3x3Definition);

%% C++ function |RefinedLeeFilter2x2| with MATLAB name |clib.gars.RefinedLeeFilter2x2|
% C++ Signature: int RefinedLeeFilter2x2(long nLooks,long height,long width,float const * c11,float const * c22,float const * c12_r,float const * c12_i,float * outC11,float * outC22,float * outC12_r,float * outC12_i)

RefinedLeeFilter2x2Definition = addFunction(libDef, ...
   "int RefinedLeeFilter2x2(long nLooks,long height,long width,float const * c11,float const * c22,float const * c12_r,float const * c12_i,float * outC11,float * outC22,float * outC12_r,float * outC12_i)", ...
   "MATLABName", "clib.gars.RefinedLeeFilter2x2", ...
   "Description", "clib.gars.RefinedLeeFilter2x2 Representation of C++ function RefinedLeeFilter2x2."); % Modify help description values as needed.
defineArgument(RefinedLeeFilter2x2Definition, "nLooks", "int32");
defineArgument(RefinedLeeFilter2x2Definition, "height", "int32");
defineArgument(RefinedLeeFilter2x2Definition, "width", "int32");
defineArgument(RefinedLeeFilter2x2Definition, "c11", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter2x2Definition, "c22", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter2x2Definition, "c12_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter2x2Definition, "c12_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter2x2Definition, "outC11", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter2x2Definition, "outC22", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter2x2Definition, "outC12_r", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(RefinedLeeFilter2x2Definition, "outC12_i", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineOutput(RefinedLeeFilter2x2Definition, "RetVal", "int32");
validate(RefinedLeeFilter2x2Definition);

%% C++ function |CloudePottier| with MATLAB name |clib.gars.CloudePottier|
% C++ Signature: void CloudePottier(long height,long width,float const * t11,float const * t22,float const * t33,float const * t12_r,float const * t13_r,float const * t23_r,float const * t12_i,float const * t13_i,float const * t23_i,float * outH,float * outAlpha,float * outA)

CloudePottierDefinition = addFunction(libDef, ...
   "void CloudePottier(long height,long width,float const * t11,float const * t22,float const * t33,float const * t12_r,float const * t13_r,float const * t23_r,float const * t12_i,float const * t13_i,float const * t23_i,float * outH,float * outAlpha,float * outA)", ...
   "MATLABName", "clib.gars.CloudePottier", ...
   "Description", "clib.gars.CloudePottier Representation of C++ function CloudePottier."); % Modify help description values as needed.
defineArgument(CloudePottierDefinition, "height", "int32");
defineArgument(CloudePottierDefinition, "width", "int32");
defineArgument(CloudePottierDefinition, "t11", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "t22", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "t33", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "t12_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "t13_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "t23_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "t12_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "t13_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "t23_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "outH", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "outAlpha", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(CloudePottierDefinition, "outA", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
validate(CloudePottierDefinition);

%% C++ function |Yamaguchi| with MATLAB name |clib.gars.Yamaguchi|
% C++ Signature: void Yamaguchi(long height,long width,float const * c11,float const * c22,float const * c33,float const * c12_r,float const * c13_r,float const * c23_r,float const * c12_i,float const * c13_i,float const * c23_i,float * outPs,float * outPd,float * outPv,float * outPh)

YamaguchiDefinition = addFunction(libDef, ...
   "void Yamaguchi(long height,long width,float const * c11,float const * c22,float const * c33,float const * c12_r,float const * c13_r,float const * c23_r,float const * c12_i,float const * c13_i,float const * c23_i,float * outPs,float * outPd,float * outPv,float * outPh)", ...
   "MATLABName", "clib.gars.Yamaguchi", ...
   "Description", "clib.gars.Yamaguchi Representation of C++ function Yamaguchi."); % Modify help description values as needed.
defineArgument(YamaguchiDefinition, "height", "int32");
defineArgument(YamaguchiDefinition, "width", "int32");
defineArgument(YamaguchiDefinition, "c11", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "c22", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "c33", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "c12_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "c13_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "c23_r", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "c12_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "c13_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "c23_i", "single", "input", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "outPs", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "outPd", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "outPv", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
defineArgument(YamaguchiDefinition, "outPh", "single", "output", ["height" "width"]); % <MLTYPE> can be "single", or "single"
validate(YamaguchiDefinition);

%% Validate the library definition
validate(libDef);

end
