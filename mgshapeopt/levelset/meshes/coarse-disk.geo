SetFactory("OpenCASCADE");
a() = ShapeFromFile("/home/wechsung/Dropbox/Documents/Uni/DPhil/ThesisNumerics/mgshapeopt/levelset/meshes/disk.step");

Mesh.CharacteristicLengthMin = 0.75;
Mesh.CharacteristicLengthMax = 0.75;
                Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Surface(4) = {1};
