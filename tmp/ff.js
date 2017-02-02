var VectorTools;
(function (VectorTools) {
    var Vector = (function () {
        function Vector(_x, _y, _z) {
            if (_x === void 0) { _x = 0; }
            if (_y === void 0) { _y = 0; }
            if (_z === void 0) { _z = 0; }
            this.x = _x;
            this.y = _y;
            this.z = _z;
        }
        Vector.prototype.dot = function (other) {
            return (this.x * other.x + this.y * other.y + this.z * other.z);
        };
        Vector.prototype.normSq = function () {
            return this.dot(this);
        };
        Vector.prototype.norm = function () {
            return Math.sqrt(this.normSq());
        };
        Vector.prototype.set = function (x, y, z) {
            this.x = x;
            this.y = y;
            this.z = z;
            return this;
        };
        Vector.sub = function (v1, v2) {
            return new Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
        };
        Vector.osub = function (u, v, out) {
            return out.set(u.x - v.x, u.y - v.y, u.z - v.z);
        };
        Vector.add = function (u, v) {
            return new Vector(u.x + v.x, u.y + v.y, u.z + v.z);
        };
        Vector.cross = function (u, v) {
            return new Vector(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
        };
        Vector.ocross = function (u, v, out) {
            return out.set(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
        };
        Vector.mult = function (v, s) {
            return new Vector(v.x * s, v.y * s, v.z * s);
        };
        Vector.omult = function (v, s, out) {
            return out.set(s * v.x, s * v.y, s * v.z);
        };
        Vector.prototype.iadd = function (v) {
            this.x += v.x;
            this.y += v.y;
            this.z += v.z;
        };
        Vector.prototype.isub = function (v) {
            this.x -= v.x;
            this.y -= v.y;
            this.z -= v.z;
        };
        Vector.prototype.imult = function (s) {
            this.x *= s;
            this.y *= s;
            this.z *= s;
        };
        Vector.prototype.izero = function () {
            this.x = this.y = this.z = 0;
        };
        return Vector;
    }());
    VectorTools.Vector = Vector;
    var VectorPool = (function () {
        function VectorPool() {
        }
        VectorPool.collect = function () {
            var v;
            while (VectorPool.requested.length > 0) {
                v = VectorPool.requested.pop();
                VectorPool.available.push(v);
            }
        };
        VectorPool.getVector = function () {
            if (VectorPool.available.length === 0) {
                VectorPool.available.push(new Vector());
            }
            var v = VectorPool.available.pop();
            VectorPool.requested.push(v);
            v.set(0, 0, 0);
            return v;
        };
        return VectorPool;
    }());
    VectorPool.available = [];
    VectorPool.requested = [];
    VectorTools.VectorPool = VectorPool;
})(VectorTools || (VectorTools = {}));
/// <reference path="VectorTools.ts" />
var Vector = VectorTools.Vector;
var VectorPool = VectorTools.VectorPool;
/// <reference path="ForceField.ts" />
var ForceField;
(function (ForceField) {
    var Angle = (function () {
        function Angle(a1, a2, a3, kf, r0) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
            this.kf = kf;
            this.r0 = r0;
        }
        Angle.eval = function (angles, coords, forces) {
            var v12 = VectorPool.getVector();
            var v32 = VectorPool.getVector();
            var f1 = VectorPool.getVector();
            var f3 = VectorPool.getVector();
            var d32Sq, d12Sq, d12xd32, fc, energy = 0.0;
            angles.forEach(function (angle) {
                Vector.osub(coords[angle.a2], coords[angle.a1], v12);
                Vector.osub(coords[angle.a2], coords[angle.a3], v32);
                d12Sq = v12.normSq();
                d32Sq = v32.normSq();
                d12xd32 = Math.sqrt(d12Sq * d32Sq);
                var ctheta = v12.dot(v32) / d12xd32;
                var dtheta = Math.acos(ctheta) - angle.r0;
                if (forces) {
                    fc = angle.kf * dtheta / Math.sqrt(1.0 - ctheta);
                    f1.set(fc * (ctheta * v12.x / d12Sq - v32.x / d12xd32), fc * (ctheta * v12.y / d12Sq - v32.y / d12xd32), fc * (ctheta * v12.z / d12Sq - v32.z / d12xd32));
                    f3.set(fc * (ctheta * v32.x / d32Sq - v12.x / d12xd32), fc * (ctheta * v32.y / d32Sq - v12.y / d12xd32), fc * (ctheta * v32.z / d32Sq - v12.z / d12xd32));
                    forces[angle.a1].iadd(f1);
                    forces[angle.a3].iadd(f3);
                    f1.iadd(f3);
                    forces[angle.a2].isub(f1);
                }
                energy += angle.kf * dtheta * dtheta;
            });
            VectorPool.collect();
            return 0.5 * energy;
        };
        return Angle;
    }());
    ForceField.Angle = Angle;
})(ForceField || (ForceField = {}));
/// <reference path="ForceField.ts" />
var ForceField;
(function (ForceField) {
    var Atom = (function () {
        function Atom(id, name, sigma, epsilon, charge) {
            this.sigma = sigma;
            this.epsilon = epsilon;
            this.charge = charge;
            this.id = id;
            this.name = name;
            this.blist = [Number.POSITIVE_INFINITY];
        }
        Atom.prototype.addBondedAtom = function (atomId) {
            if (this.blist.indexOf(atomId) < 0) {
                this.blist.push(atomId);
            }
            this.blist.sort(function (a, b) { return a - b; });
        };
        Atom.eval = function (atoms, coords, forces) {
            var vij = VectorPool.getVector();
            var nAtoms = atoms.length;
            var bidx, bid;
            var atomi, atomj;
            var dijSq, sij, eij, fc, lj6, Ecoul, Evdw, EnergyVDW = 0.0, EnergyCoul = 0.0;
            for (var i = 0; i < nAtoms - 1; i++) {
                bidx = 0;
                atomi = atoms[i];
                bid = atomi[bidx];
                for (var j = i + 1; j < nAtoms; j++) {
                    if (j === bid) {
                        bid = atomi[++bidx];
                        continue;
                    }
                    Vector.osub(coords[j], coords[i], vij);
                    dijSq = vij.normSq();
                    atomj = atoms[j];
                    sij = atomi.sigma + atomj.sigma;
                    eij = atomi.epsilon * atomj.epsilon;
                    lj6 = Math.pow(sij * sij / dijSq, 3);
                    Evdw = 4.0 * eij * (lj6 * lj6 - lj6);
                    Ecoul = atomi.charge * atomj.charge / Math.sqrt(dijSq);
                    if (forces) {
                        fc = -(24 * eij * (2 * lj6 * lj6 - lj6) + Ecoul) / dijSq;
                        vij.imult(fc);
                        forces[i].iadd(vij);
                        forces[j].isub(vij);
                    }
                    EnergyCoul += Ecoul;
                    EnergyVDW += Evdw;
                }
            }
            VectorPool.collect();
            return { coulomb: EnergyCoul, vdw: EnergyVDW };
        };
        return Atom;
    }());
    ForceField.Atom = Atom;
})(ForceField || (ForceField = {}));
/// <reference path="ForceField.ts" />
var ForceField;
(function (ForceField) {
    var Bond = (function () {
        function Bond(a1, a2, kf, r0) {
            this.a1 = a1;
            this.a2 = a2;
            this.kf = kf;
            this.r0 = r0;
        }
        Bond.eval = function (bonds, coords, forces) {
            var v12 = VectorPool.getVector(); // new Vector();
            var d12, dr, energy = 0.0;
            bonds.forEach(function (bond) {
                Vector.osub(coords[bond.a2], coords[bond.a1], v12);
                d12 = v12.norm();
                dr = d12 - bond.r0;
                if (forces) {
                    v12.imult(bond.kf * dr / d12);
                    forces[bond.a1].iadd(v12);
                    forces[bond.a2].isub(v12);
                }
                energy += bond.kf * dr * dr;
            });
            VectorPool.collect();
            return 0.5 * energy;
        };
        return Bond;
    }());
    ForceField.Bond = Bond;
})(ForceField || (ForceField = {}));
/// <reference path="ForceField.ts" />
var ForceField;
(function (ForceField) {
    var DihedralCosine = (function () {
        function DihedralCosine(a1, a2, a3, a4, kf, pn, ph) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
            this.a4 = a4;
            this.kf = kf;
            this.pn = pn;
            this.ph = ph;
        }
        DihedralCosine.eval = function (dihedrals, coords, forces) {
            var v12 = VectorPool.getVector();
            var v32 = VectorPool.getVector();
            var v34 = VectorPool.getVector();
            var tmp = VectorPool.getVector();
            var M = VectorPool.getVector();
            var N = VectorPool.getVector();
            var f1 = VectorPool.getVector();
            var f2 = VectorPool.getVector();
            var f3 = VectorPool.getVector();
            var f4 = VectorPool.getVector();
            var d32Sq, d32, phi, dVdphi_x_d32, energy = 0.0;
            dihedrals.forEach(function (dihd) {
                Vector.osub(coords[dihd.a1], coords[dihd.a2], v12);
                Vector.osub(coords[dihd.a3], coords[dihd.a2], v32);
                Vector.osub(coords[dihd.a3], coords[dihd.a4], v34);
                Vector.ocross(v12, v32, M);
                Vector.ocross(v32, v34, N);
                d32Sq = v32.normSq();
                d32 = Math.sqrt(d32Sq);
                phi = Math.atan2(d32 * v12.dot(N), M.dot(N)) - Math.PI;
                if (forces) {
                    dVdphi_x_d32 = dihd.kf * dihd.pn * Math.sin(dihd.pn * phi) * d32;
                    Vector.omult(M, -dVdphi_x_d32 / M.dot(M), f1); // F1
                    Vector.omult(N, dVdphi_x_d32 / N.dot(N), f4); // F4
                    f3.set(-f4.x, -f4.y, -f4.z);
                    f3.isub(Vector.omult(f1, v12.dot(v32) / d32Sq, tmp));
                    f3.iadd(Vector.omult(f4, v34.dot(v32) / d32Sq, tmp));
                    tmp.set(-f1.x - f3.x - f4.x, -f1.y - f3.y - f4.y, -f1.z - f3.z - f4.z);
                    forces[dihd.a1].iadd(f1);
                    forces[dihd.a2].iadd(tmp);
                    forces[dihd.a3].iadd(f3);
                    forces[dihd.a4].iadd(f4);
                }
                energy += dihd.kf * (1.0 + Math.cos(dihd.pn * phi - dihd.ph));
            });
            VectorPool.collect();
            return energy;
        };
        return DihedralCosine;
    }());
    ForceField.DihedralCosine = DihedralCosine;
})(ForceField || (ForceField = {}));
/// <reference path="ForceField.ts" />
var ForceField;
(function (ForceField) {
    var DihedralRB = (function () {
        function DihedralRB(a1, a2, a3, a4, c0, c1, c2, c3, c4, c5) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
            this.a4 = a4;
            this.c0 = c0;
            this.c1 = c1;
            this.c2 = c2;
            this.c3 = c3;
            this.c4 = c4;
            this.c5 = c5;
        }
        DihedralRB.eval = function (dihedrals, coords, forces) {
            var v12 = VectorPool.getVector();
            var v32 = VectorPool.getVector();
            var v34 = VectorPool.getVector();
            var tmp = VectorPool.getVector();
            var f1 = VectorPool.getVector();
            var f2 = VectorPool.getVector();
            var f3 = VectorPool.getVector();
            var f4 = VectorPool.getVector();
            var M = VectorPool.getVector();
            var N = VectorPool.getVector();
            var d32Sq, d32, phi, cpsi, cphi, dVdphi_x_d32, energy = 0.0;
            dihedrals.forEach(function (dihd) {
                Vector.osub(coords[dihd.a1], coords[dihd.a2], v12);
                Vector.osub(coords[dihd.a3], coords[dihd.a2], v32);
                Vector.osub(coords[dihd.a3], coords[dihd.a4], v34);
                Vector.ocross(v12, v32, M);
                Vector.ocross(v32, v34, N);
                d32Sq = v32.normSq();
                d32 = Math.sqrt(d32Sq);
                phi = Math.atan2(d32 * v12.dot(N), M.dot(N));
                cpsi = Math.cos(phi - Math.PI);
                energy += dihd.c0 + cpsi * (dihd.c1 + cpsi * (dihd.c2 + cpsi * (dihd.c3 + cpsi * (dihd.c4 + cpsi * dihd.c5))));
                if (forces) {
                    cphi = -cpsi; // cos(phi-PI) = cos(psi) = -cos(phi)
                    dVdphi_x_d32 = Math.sin(phi) * d32 * (dihd.c1 + cphi * (-2.0 * dihd.c2 + cphi * (3.0 * dihd.c3 + cphi * (-4.0 * dihd.c4 + cphi * 5.0 * dihd.c5))));
                    Vector.omult(M, -dVdphi_x_d32 / M.dot(M), f1); // F1
                    Vector.omult(N, dVdphi_x_d32 / N.dot(N), f4); // F4
                    f3.set(-f4.x, -f4.y, -f4.z);
                    f3.isub(Vector.omult(f1, v12.dot(v32) / d32Sq, tmp));
                    f3.iadd(Vector.omult(f4, v34.dot(v32) / d32Sq, tmp));
                    tmp.set(-f1.x - f3.x - f4.x, -f1.y - f3.y - f4.y, -f1.z - f3.z - f4.z);
                    forces[dihd.a1].iadd(f1);
                    forces[dihd.a2].iadd(tmp);
                    forces[dihd.a3].iadd(f3);
                    forces[dihd.a4].iadd(f4);
                }
            });
            VectorPool.collect();
            return energy;
        };
        return DihedralRB;
    }());
    ForceField.DihedralRB = DihedralRB;
})(ForceField || (ForceField = {}));
/// <reference path="ForceField.ts" />
var MolecularMechanics;
(function (MolecularMechanics) {
    var System = (function () {
        function System() {
            this.dihedralsCos = [];
            this.dihedralsRB = [];
            this.angles = [];
            this.bonds = [];
            this.atoms = [];
        }
        System.prototype.addBond = function (a1, a2, kf, r0) {
            var item = new ForceField.Bond(a1, a2, kf, r0);
            this.bonds.push(item);
            return item;
        };
        System.prototype.addAngle = function (a1, a2, a3, kf, r0) {
            var item = new ForceField.Angle(a1, a2, a3, kf, r0);
            this.angles.push(item);
            return item;
        };
        System.prototype.addCosineDihedral = function (a1, a2, a3, a4, kf, pn, ph) {
            var item = new ForceField.DihedralCosine(a1, a2, a3, a4, kf, pn, ph);
            this.dihedralsCos.push(item);
            return item;
        };
        System.prototype.addRBDihedral = function (a1, a2, a3, a4, c0, c1, c2, c3, c4, c5) {
            var item = new ForceField.DihedralRB(a1, a2, a3, a4, c0, c1, c2, c3, c4, c5);
            this.dihedralsRB.push(item);
            return item;
        };
        System.prototype.addAtom = function (id, name, sigma, epsilon, charge) {
            var item = new ForceField.Atom(id, name, sigma, epsilon, charge);
            this.atoms.push(item);
            return item;
        };
        System.prototype.eval = function (coords, forces) {
            var Ebonded = 0;
            var Enonbonded;
            Ebonded += ForceField.Bond.eval(this.bonds, coords, forces);
            Ebonded += ForceField.Angle.eval(this.angles, coords, forces);
            Ebonded += ForceField.DihedralRB.eval(this.dihedralsRB, coords, forces);
            Ebonded += ForceField.DihedralCosine.eval(this.dihedralsCos, coords, forces);
            Enonbonded = ForceField.Atom.eval(this.atoms, coords, forces);
            return Ebonded + Enonbonded.coulomb + Enonbonded.vdw;
        };
        return System;
    }());
    MolecularMechanics.System = System;
})(MolecularMechanics || (MolecularMechanics = {}));
/// <reference path="System.ts"/>
var fs = require('fs');
var ffdata = JSON.parse(fs.readFileSync('ffdat.js', 'utf8'));
var mm = new MolecularMechanics.System();
function deg2rad(x) {
    return x * Math.PI / 180;
}
ffdata.atoms.forEach(function (it) {
    var item = mm.addAtom(it.id, it.name, it.sigma, it.epsilon, it.charge);
    console.log(item);
});
ffdata.bonds.forEach(function (it) {
    var item = mm.addBond(it.a1, it.a2, it.kf, it.r0);
    console.log(item);
});
ffdata.angles.forEach(function (it) {
    var item = mm.addAngle(it.a1, it.a2, it.a3, it.kf, deg2rad(it.r0));
    console.log(item);
});
ffdata.dihedralsRB.forEach(function (it) {
    var item = mm.addRBDihedral(it.a1, it.a2, it.a3, it.a4, it.c0, it.c1, it.c2, it.c3, it.c4, it.c5);
    console.log(item);
});
ffdata.dihedralsCos.forEach(function (it) {
    var item = mm.addCosineDihedral(it.a1, it.a2, it.a3, it.a4, it.kf, it.pn, deg2rad(it.ph));
    console.log(item);
});
