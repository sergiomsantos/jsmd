/// <reference path="System.ts"/>

var fs = require('fs');
var ffdata = JSON.parse(fs.readFileSync('ffdat.js', 'utf8'));

var mm = new MolecularMechanics.System();

function deg2rad(x: number): number {
  return x*Math.PI/180;
}

ffdata.atoms.forEach(it => {
  var item = mm.addAtom(it.id, it.name, it.sigma, it.epsilon, it.charge);
  console.log(item);
});

ffdata.bonds.forEach(it => {
  var item = mm.addBond(it.a1, it.a2, it.kf, it.r0);
  console.log(item);
});

ffdata.angles.forEach(it => {
  var item = mm.addAngle(it.a1, it.a2, it.a3, it.kf, deg2rad(it.r0));
  console.log(item);
});

ffdata.dihedralsRB.forEach(it => {
  var item = mm.addRBDihedral(it.a1, it.a2, it.a3, it.a4,
                  it.c0, it.c1, it.c2, it.c3, it.c4, it.c5);
  console.log(item);
});

ffdata.dihedralsCos.forEach(it => {
  var item = mm.addCosineDihedral(it.a1, it.a2, it.a3, it.a4, it.kf, it.pn, deg2rad(it.ph));
  console.log(item);
});
