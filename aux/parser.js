var fs = require('fs');
var obj = JSON.parse(fs.readFileSync('rawff.json', 'utf8'));

var top = obj.Topology.Molecule;

var parsed = {};
parsed.atoms = [];
top.Nonbonded.atom.forEach(function(item) {
  parsed.atoms.push({
    sigma: parseFloat(item.ljs),
    epsilon: parseFloat(item.lje),
    charge: parseFloat(item.q),
    name: item.name,
    id: parseInt(item.ID)
  });
});


parsed.bonds = [];
top.Bonds.bond.forEach(function(item) {
  parsed.bonds.push({
    a1: parseInt(item.a1),
    a2: parseInt(item.a2),
    kf: parseFloat(item.kf),
    r0: parseFloat(item.x0)
  });
});

parsed.angles = [];
top.Angles.angle.forEach(function(item) {
  parsed.angles.push({
    a1: parseInt(item.a1),
    a2: parseInt(item.a2),
    a3: parseInt(item.a3),
    kf: parseFloat(item.kf),
    r0: parseFloat(item.x0)
  });
});


parsed.dihedralsRB  = [];
parsed.dihedralsCos = [];
top.Dihedrals.dihedral.forEach(function(item) {
  if (item.type==='cos') {
    parsed.dihedralsCos.push({
      a1: parseInt(item.a1),
      a2: parseInt(item.a2),
      a3: parseInt(item.a3),
      a4: parseInt(item.a4),
      kf: parseFloat(item.kf),
      ph: parseFloat(item.ph),
      pn: parseInt(item.pn)
    });
  } else {
    parsed.dihedralsRB.push({
      a1: parseInt(item.a1),
      a2: parseInt(item.a2),
      a3: parseInt(item.a3),
      a4: parseInt(item.a4),
      c0: parseFloat(item.c0)||0,
      c1: parseFloat(item.c1)||0,
      c2: parseFloat(item.c2)||0,
      c3: parseFloat(item.c3)||0,
      c4: parseFloat(item.c4)||0,
      c5: parseFloat(item.c5)||0,
    });
  }
});

console.log(JSON.stringify(parsed));
// obj = JSON.parse(fs.readFileSync('ffdat.js', 'utf8'));