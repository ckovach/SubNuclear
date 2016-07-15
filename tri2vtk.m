
function vtk = tri2vtk(Tr)

% Convert TriRep to vtk struct

vtk.format = 'float';

vtk.tri = Tr.Triangulation;
vtk.n = length(Tr.X);
vtk.vert = Tr.X;
