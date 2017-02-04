var __extends = (this && this.__extends) || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
};
var Vectors;
(function (Vectors) {
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
        Vector.prototype.copy = function (v) {
            this.x = v.x;
            this.y = v.y;
            this.z = v.z;
            return this;
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
        Vector.oadd = function (u, v, out) {
            return out.set(u.x + v.x, u.y + v.y, u.z + v.z);
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
            return this;
        };
        Vector.prototype.isub = function (v) {
            this.x -= v.x;
            this.y -= v.y;
            this.z -= v.z;
            return this;
        };
        Vector.prototype.imult = function (s) {
            this.x *= s;
            this.y *= s;
            this.z *= s;
            return this;
        };
        Vector.prototype.izero = function () {
            this.x = this.y = this.z = 0;
            return this;
        };
        Vector.prototype.free = function () {
            VectorPool.uptake(this);
        };
        Vector.getVector = function (n) {
            if (n) {
                var vs = [];
                while ((n--) > 0)
                    vs.push(new Vector());
                return vs;
            }
            return (new Vector());
        };
        return Vector;
    }());
    Vectors.Vector = Vector;
    var VectorPool = (function () {
        function VectorPool() {
        }
        VectorPool.collect = function () {
            var v;
            while (VectorPool.requested.length > 0) {
                v = VectorPool.requested.pop();
                VectorPool.allocated.push(v);
            }
        };
        VectorPool.uptake = function (v) {
            this.allocated.push(v);
        };
        VectorPool.getVector = function (n) {
            if (n) {
                var vs = [];
                while ((n--) > 0)
                    vs.push(VectorPool.getVector());
                return vs;
            }
            else {
                if (VectorPool.allocated.length === 0) {
                    VectorPool.allocated.push(new Vector());
                    VectorPool.counter++;
                }
                var v = VectorPool.allocated.pop();
                VectorPool.requested.push(v);
                v.izero();
                return v;
            }
        };
        return VectorPool;
    }());
    VectorPool.counter = 0;
    VectorPool.allocated = [];
    VectorPool.requested = [];
    Vectors.VectorPool = VectorPool;
})(Vectors || (Vectors = {}));
var Vector = Vectors.Vector;
var VectorPool = Vectors.VectorPool;
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
            v12.free();
            v32.free();
            f1.free();
            f3.free();
            return 0.5 * energy;
        };
        return Angle;
    }());
    ForceField.Angle = Angle;
})(ForceField || (ForceField = {}));
var ForceField;
(function (ForceField) {
    var Atom = (function () {
        function Atom(id, name, sigma, epsilon, charge) {
            this.sigma = sigma;
            this.epsilon = epsilon;
            this.charge = charge;
            this.id = id - 1;
            this.name = name;
            this.blist = [Number.POSITIVE_INFINITY];
            this.mass = Atom.MASS[name.replace(/[0-9]+/, '')];
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
                bidx = -1;
                atomi = atoms[i];
                while (atomi.blist[++bidx] < i) { }
                bid = atomi.blist[bidx];
                for (var j = i + 1; j < nAtoms; j++) {
                    if (j === bid) {
                        bid = atomi.blist[++bidx];
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
                        fc = (24 * eij * (2 * lj6 * lj6 - lj6) - Ecoul) / dijSq;
                        vij.imult(-fc);
                        forces[i].iadd(vij);
                        forces[j].isub(vij);
                    }
                    EnergyCoul += Ecoul;
                    EnergyVDW += Evdw;
                }
            }
            vij.free();
            return { coulomb: EnergyCoul, vdw: EnergyVDW };
        };
        return Atom;
    }());
    Atom.MASS = {
        'H': 1,
        'C': 12,
        'N': 14,
        'O': 16
    };
    ForceField.Atom = Atom;
})(ForceField || (ForceField = {}));
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
            var v12 = VectorPool.getVector();
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
            v12.free();
            return 0.5 * energy;
        };
        return Bond;
    }());
    ForceField.Bond = Bond;
})(ForceField || (ForceField = {}));
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
                    Vector.omult(M, -dVdphi_x_d32 / M.dot(M), f1);
                    Vector.omult(N, dVdphi_x_d32 / N.dot(N), f4);
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
            v12.free();
            v32.free();
            v34.free();
            tmp.free();
            M.free();
            N.free();
            f1.free();
            f2.free();
            f3.free();
            f4.free();
            return energy;
        };
        return DihedralCosine;
    }());
    ForceField.DihedralCosine = DihedralCosine;
})(ForceField || (ForceField = {}));
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
                    cphi = -cpsi;
                    dVdphi_x_d32 = Math.sin(phi) * d32 * (dihd.c1 + cphi * (-2.0 * dihd.c2 + cphi * (3.0 * dihd.c3 + cphi * (-4.0 * dihd.c4 + cphi * 5.0 * dihd.c5))));
                    Vector.omult(M, -dVdphi_x_d32 / M.dot(M), f1);
                    Vector.omult(N, dVdphi_x_d32 / N.dot(N), f4);
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
            v12.free();
            v32.free();
            v34.free();
            tmp.free();
            f1.free();
            f2.free();
            f3.free();
            f4.free();
            M.free();
            N.free();
            return energy;
        };
        return DihedralRB;
    }());
    ForceField.DihedralRB = DihedralRB;
})(ForceField || (ForceField = {}));
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
        System.prototype.setCoordinates = function (coords) {
            this.x = coords;
        };
        System.prototype.getCoordinates = function () {
            return this.x;
        };
        System.prototype.setForces = function (forces) {
            this.f = forces;
        };
        System.prototype.getForces = function () {
            return this.f;
        };
        System.prototype.setVelocities = function (v) {
            this.v = v;
        };
        System.prototype.getVelocities = function () {
            return this.v;
        };
        System.prototype.getAtoms = function () {
            return this.atoms;
        };
        System.prototype.countDegreesOfFreedom = function () {
            return 3 * this.atoms.length - 6;
        };
        System.prototype.loadState = function (state) {
            var Natoms = this.getAtomCount();
            console.log(Natoms);
            this.x = Vector.getVector(Natoms);
            this.f = Vector.getVector(Natoms);
            this.v = Vector.getVector(Natoms);
            console.log(this.x);
            console.log(this.f);
            console.log(this.v);
            if (state.coords) {
                this.x.forEach(function (x, i) {
                    x.set(state.coords[i][0], state.coords[i][1], state.coords[i][2]);
                });
            }
            if (state.forces) {
                this.f.forEach(function (f, i) {
                    f.set(state.forces[i][0], state.forces[i][1], state.forces[i][2]);
                });
            }
            if (state.velocities) {
                this.v.forEach(function (v, i) {
                    v.set(state.velocities[i][0], state.velocities[i][1], state.velocities[i][2]);
                });
            }
        };
        System.prototype.loadParamsFromJSON = function (obj) {
            var _this = this;
            var deg2rad = Math.PI / 180.0;
            if (obj.atoms)
                obj.atoms.forEach(function (it) {
                    console.log(_this.addAtom(it.id, it.name, it.sigma, it.epsilon, it.charge));
                });
            if (obj.bonds)
                obj.bonds.forEach(function (it) {
                    console.log(_this.addBond(it.a1, it.a2, it.kf, it.r0));
                });
            if (obj.angles)
                obj.angles.forEach(function (it) {
                    console.log(_this.addAngle(it.a1, it.a2, it.a3, it.kf, deg2rad * it.r0));
                });
            if (obj.dihedralsRB)
                obj.dihedralsRB.forEach(function (it) {
                    console.log(_this.addRBDihedral(it.a1, it.a2, it.a3, it.a4, it.c0, it.c1, it.c2, it.c3, it.c4, it.c5));
                });
            if (obj.dihedralsCos)
                obj.dihedralsCos.forEach(function (it) {
                    console.log(_this.addCosineDihedral(it.a1, it.a2, it.a3, it.a4, it.kf, it.pn, deg2rad * it.ph));
                });
        };
        System.prototype.getAtomCount = function () {
            return this.atoms.length;
        };
        System.prototype.addBond = function (a1, a2, kf, r0) {
            var item = new ForceField.Bond(a1, a2, kf, r0);
            this.bonds.push(item);
            this.atoms[a1].addBondedAtom(a2);
            this.atoms[a2].addBondedAtom(a1);
            return item;
        };
        System.prototype.addAngle = function (a1, a2, a3, kf, r0) {
            var item = new ForceField.Angle(a1, a2, a3, kf, r0);
            this.angles.push(item);
            this.atoms[a1].addBondedAtom(a2);
            this.atoms[a1].addBondedAtom(a3);
            this.atoms[a2].addBondedAtom(a1);
            this.atoms[a2].addBondedAtom(a3);
            this.atoms[a3].addBondedAtom(a1);
            this.atoms[a3].addBondedAtom(a2);
            return item;
        };
        System.prototype.tmp = function (i, j, k, l) {
            this.atoms[i].addBondedAtom(j);
            this.atoms[i].addBondedAtom(k);
            this.atoms[i].addBondedAtom(l);
            this.atoms[j].addBondedAtom(i);
            this.atoms[j].addBondedAtom(k);
            this.atoms[j].addBondedAtom(l);
            this.atoms[k].addBondedAtom(i);
            this.atoms[k].addBondedAtom(j);
            this.atoms[k].addBondedAtom(l);
            this.atoms[l].addBondedAtom(i);
            this.atoms[l].addBondedAtom(j);
            this.atoms[l].addBondedAtom(k);
        };
        System.prototype.addCosineDihedral = function (a1, a2, a3, a4, kf, pn, ph) {
            var item = new ForceField.DihedralCosine(a1, a2, a3, a4, kf, pn, ph);
            this.dihedralsCos.push(item);
            this.tmp(a1, a2, a3, a4);
            return item;
        };
        System.prototype.addRBDihedral = function (a1, a2, a3, a4, c0, c1, c2, c3, c4, c5) {
            var item = new ForceField.DihedralRB(a1, a2, a3, a4, c0, c1, c2, c3, c4, c5);
            this.dihedralsRB.push(item);
            this.tmp(a1, a2, a3, a4);
            return item;
        };
        System.prototype.addAtom = function (id, name, sigma, epsilon, charge) {
            var sqrtfSI = 37.27405062506623;
            var item = new ForceField.Atom(id, name, 0.5 * sigma, Math.sqrt(epsilon), sqrtfSI * charge);
            this.atoms.push(item);
            return item;
        };
        System.prototype.eval = function (coords, forces) {
            var Ebonded = 0;
            var Enonbonded;
            if (forces)
                forces.forEach(function (force) { return force.izero(); });
            if (this.bonds.length > 0)
                Ebonded += ForceField.Bond.eval(this.bonds, coords, forces);
            if (this.angles.length > 0)
                Ebonded += ForceField.Angle.eval(this.angles, coords, forces);
            if (this.dihedralsRB.length > 0)
                Ebonded += ForceField.DihedralRB.eval(this.dihedralsRB, coords, forces);
            if (this.dihedralsCos.length > 0)
                Ebonded += ForceField.DihedralCosine.eval(this.dihedralsCos, coords, forces);
            if (this.atoms.length > 0)
                Enonbonded = ForceField.Atom.eval(this.atoms, coords, forces);
            return Ebonded + Enonbonded.coulomb + Enonbonded.vdw;
        };
        return System;
    }());
    MolecularMechanics.System = System;
})(MolecularMechanics || (MolecularMechanics = {}));
var boltzman = 0.0083145112119;
var Thermostat;
(function (Thermostat) {
    var BaseThermostat = (function () {
        function BaseThermostat(temperature) {
            this.temperature = temperature;
        }
        BaseThermostat.prototype.assignVelocities = function (system) {
            var randn = function () {
                var u = 1 - Math.random();
                var v = 1 - Math.random();
                return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
            };
            var masses = system.getAtoms().map(function (atom) { return atom.mass; });
            var velocities = system.getVelocities();
            var f = this.temperature * boltzman;
            velocities.forEach(function (vel, i) { return vel.set(f * masses[i] * randn(), f * masses[i] * randn(), f * masses[i] * randn()); });
            var Ti = this.getInstantTemperature(velocities, masses, system.countDegreesOfFreedom());
            f = Math.sqrt(this.temperature / Ti);
            velocities.forEach(function (vel) { return vel.imult(f); });
            var vcom = VectorPool.getVector();
            var vtmp = VectorPool.getVector();
            vcom.izero();
            velocities.forEach(function (vel, i) {
                Vector.omult(vel, masses[i], vtmp);
                vcom.iadd(vtmp);
            });
            var totalMass = masses.reduce(function (a, b) { return a + b; }, 0);
            vcom.imult(1.0 / totalMass);
            velocities.forEach(function (v) { return v.isub(vcom); });
        };
        BaseThermostat.prototype.getInstantTemperature = function (velocities, masses, Ndf) {
            var Ti = masses.reduce(function (ekin, mass, i) { return ekin + mass * velocities[i].normSq(); }, 0.0);
            return Ti / (Ndf * boltzman);
        };
        return BaseThermostat;
    }());
    var Berendsen = (function (_super) {
        __extends(Berendsen, _super);
        function Berendsen(params) {
            var _this = _super.call(this, params.temperature) || this;
            _this.tcoupl = params.tcoupl;
            _this.dt = params.dt;
            return _this;
        }
        Berendsen.prototype.update = function (velocities, instantTemperature) {
            var scale = Math.sqrt(1.0 + this.dt * (this.temperature / instantTemperature - 1.0) / this.tcoupl);
            velocities.forEach(function (v) { return v.imult(scale); });
        };
        return Berendsen;
    }(BaseThermostat));
    Thermostat.Berendsen = Berendsen;
    var VRescale = (function (_super) {
        __extends(VRescale, _super);
        function VRescale(params) {
            return _super.call(this, params.temperature) || this;
        }
        VRescale.prototype.update = function (velocities, instantTemperature) {
            var scale = Math.sqrt(this.temperature / instantTemperature);
            velocities.forEach(function (v) { return v.imult(scale); });
        };
        return VRescale;
    }(BaseThermostat));
    Thermostat.VRescale = VRescale;
})(Thermostat || (Thermostat = {}));
var Driver;
(function (Driver) {
    var LeapfrogIntegrator = (function () {
        function LeapfrogIntegrator(params) {
            this.nsteps = params.nsteps;
            this.timestep = params.dt;
            this.momRemover = new MomentumRemover({ removeLinear: true, removeAngular: true, freq: 10 });
        }
        LeapfrogIntegrator.prototype.setCallback = function (fcn, callfreq) {
            this.callback = fcn;
            this.callbackfreq = callfreq;
        };
        LeapfrogIntegrator.prototype.setThermostat = function (thermostat) {
            this.thermostat = thermostat;
        };
        LeapfrogIntegrator.prototype.run = function (system) {
            var _this = this;
            var masses = system.getAtoms().map(function (atom) { return atom.mass; });
            var totalMass = masses.reduce(function (a, b) { return a + b; }, 0);
            var Ndf = system.countDegreesOfFreedom();
            var velocities = system.getVelocities();
            var coords = system.getCoordinates();
            var Natoms = system.getAtomCount();
            var vtmp = VectorPool.getVector();
            var forces = system.getForces();
            var instantTemp;
            var halfdtOverMass = masses.map(function (m) { return 0.5 * _this.timestep / m; });
            var forces_t = VectorPool.getVector(Natoms);
            var energy = system.eval(coords, forces);
            instantTemp = this.thermostat.getInstantTemperature(velocities, masses, Ndf);
            var step = -1;
            while ((++step) < this.nsteps) {
                this.thermostat.update(velocities, instantTemp);
                for (var i = 0; i < Natoms; i++) {
                    Vector.omult(forces[i], halfdtOverMass[i], vtmp);
                    vtmp.iadd(velocities[i]);
                    vtmp.imult(this.timestep);
                    coords[i].iadd(vtmp);
                    forces_t[i].copy(forces[i]);
                }
                energy = system.eval(coords, forces);
                for (var i = 0; i < Natoms; i++) {
                    Vector.oadd(forces_t[i], forces[i], vtmp);
                    vtmp.imult(halfdtOverMass[i]);
                    velocities[i].iadd(vtmp);
                }
                instantTemp = this.thermostat.getInstantTemperature(velocities, masses, Ndf);
                if (this.momRemover && (step % 20 == 0)) {
                    this.momRemover.remove(coords, velocities, masses);
                }
            }
            if (this.callback) {
                this.callback(system, step);
            }
            vtmp.free();
            while (forces_t.length > 0)
                forces_t.pop().free();
            return -1;
        };
        return LeapfrogIntegrator;
    }());
    Driver.LeapfrogIntegrator = LeapfrogIntegrator;
    var MomentumRemover = (function () {
        function MomentumRemover(prm) {
            this.removeAngular = prm.removeAngular;
            this.removeLinear = prm.removeLinear;
            this.freq = prm.freq;
        }
        MomentumRemover.prototype.remove = function (xyz, vel, masses) {
            var p0 = VectorPool.getVector();
            var j0 = VectorPool.getVector();
            var gj = VectorPool.getVector();
            var gx = VectorPool.getVector();
            var x0 = VectorPool.getVector();
            var gv = VectorPool.getVector();
            var gp = VectorPool.getVector();
            var gw = VectorPool.getVector();
            var dx = VectorPool.getVector();
            var dv = VectorPool.getVector();
            var gi = VectorPool.getVector(3);
            var Icm = VectorPool.getVector(3);
            var nAtoms = xyz.length;
            var m0;
            var M = 0.0;
            for (var i = 0; i < nAtoms; i++) {
                m0 = masses[i];
                M += m0;
                Vector.omult(vel[i], m0, p0);
                gp.iadd(p0);
                Vector.ocross(xyz[i], vel[i], j0);
                j0.imult(m0);
                gj.iadd(j0);
                gx.set(gx.x + m0 * xyz[i].x, gx.y + m0 * xyz[i].y, gx.z + m0 * xyz[i].z);
                this.update_tensor(xyz[i], m0, gi);
            }
            Vector.omult(gp, 1.0 / M, gv);
            gx.imult(1.0 / M);
            var jcm = VectorPool.getVector();
            Vector.ocross(gx, gv, jcm);
            jcm.imult(M);
            gj.isub(jcm);
            this.update_tensor(gx, M, Icm);
            gi[0].isub(Icm[0]);
            gi[1].isub(Icm[1]);
            gi[2].isub(Icm[2]);
            this.get_minv(gi, Icm);
            gw.set(Icm[0].dot(gj), Icm[1].dot(gj), Icm[2].dot(gj));
            for (var i = 0; i < nAtoms; i++) {
                vel[i].isub(gv);
                Vector.osub(xyz[i], gx, dx);
                Vector.ocross(gw, dx, dv);
                vel[i].isub(dv);
            }
            p0.free();
            j0.free();
            gj.free();
            gx.free();
            x0.free();
            gv.free();
            gp.free();
            gw.free();
            dx.free();
            dv.free();
            while (gi.length > 0)
                gi.pop().free();
            while (Icm.length > 0)
                Icm.pop().free();
            jcm.free();
        };
        MomentumRemover.prototype.update_tensor = function (vec, s, tensor) {
            var xy = vec.x * vec.y * s;
            var xz = vec.x * vec.z * s;
            var yz = vec.y * vec.z * s;
            tensor[0].x += vec.x * vec.x * s;
            tensor[1].y += vec.y * vec.y * s;
            tensor[2].z += vec.z * vec.z * s;
            tensor[0].y += xy;
            tensor[1].x += xy;
            tensor[0].z += xz;
            tensor[2].x += xz;
            tensor[1].z += yz;
            tensor[2].y += yz;
        };
        MomentumRemover.prototype.get_minv = function (A, B) {
            var tmp = VectorPool.getVector(3);
            tmp[0].x = A[1].y + A[2].z;
            tmp[1].x = -A[0].y;
            tmp[2].x = -A[0].z;
            tmp[0].y = -A[0].y;
            tmp[1].y = A[0].x + A[2].z;
            tmp[2].y = -A[1].z;
            tmp[0].z = -A[0].z;
            tmp[1].z = -A[1].z;
            tmp[2].z = A[0].x + A[1].y;
            var rfac = (tmp[0].x + tmp[1].y + tmp[2].z) / 3;
            var fac = 1 / rfac;
            tmp[0].imult(fac);
            tmp[1].imult(fac);
            tmp[2].imult(fac);
            this.m_inv(tmp[0].x, tmp[0].y, tmp[0].z, tmp[1].x, tmp[1].y, tmp[1].z, tmp[2].x, tmp[2].y, tmp[2].z, B);
            B[0].imult(fac);
            B[1].imult(fac);
            B[2].imult(fac);
            while (tmp.length > 0)
                tmp.pop().free();
        };
        MomentumRemover.prototype.m_inv = function (a11, a12, a13, a21, a22, a23, a31, a32, a33, inv) {
            var det = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) + a31 * (a23 * a12 - a22 * a13);
            inv[0].set(a33 * a22 - a32 * a23, -(a33 * a12 - a32 * a13), a23 * a12 - a22 * a13);
            inv[1].set(-(a33 * a21 - a31 * a23), a33 * a11 - a31 * a13, -(a23 * a11 - a21 * a13));
            inv[2].set(a32 * a21 - a31 * a22, -(a32 * a11 - a31 * a12), a22 * a11 - a21 * a12);
            inv[0].imult(1 / det);
            inv[1].imult(1 / det);
            inv[2].imult(1 / det);
        };
        return MomentumRemover;
    }());
    Driver.MomentumRemover = MomentumRemover;
    var SteepestDescent = (function () {
        function SteepestDescent() {
            this.nSteps = 10000;
            this.Etol = 0.00001;
            this.stepsize = 1;
        }
        SteepestDescent.prototype.setCallback = function (fcn) {
            this.callback = fcn;
        };
        SteepestDescent.prototype.run = function (system) {
            var Eold;
            var Enew = +Infinity;
            var forces = system.getForces();
            var coords = system.getCoordinates();
            var normSqr;
            var stepsize;
            normSqr = forces.reduce(function (old, cur, i) {
                var nSqr = cur.normSq();
                return nSqr > old ? nSqr : old;
            }, -1);
            var step = 0;
            while ((++step) < this.nSteps) {
                Eold = Enew;
                Enew = system.eval(coords, forces);
                console.log(step, 'stepsize', this.stepsize, Enew - Eold);
                if (Math.abs(Eold - Enew) < this.Etol)
                    break;
                if (Enew > Eold)
                    this.stepsize *= 0.5;
                else
                    this.stepsize *= 1.05;
                normSqr = forces.reduce(function (old, cur) {
                    var nSqr = cur.normSq();
                    return (nSqr > old ? nSqr : old);
                }, -1);
                stepsize = this.stepsize / Math.sqrt(normSqr);
                for (var i = 0; i < coords.length; i++) {
                    forces[i].imult(stepsize);
                    coords[i].iadd(forces[i]);
                }
                if (this.callback)
                    this.callback(system, Enew);
            }
            return Enew;
        };
        return SteepestDescent;
    }());
    Driver.SteepestDescent = SteepestDescent;
})(Driver || (Driver = {}));
//# sourceMappingURL=ff.js.map