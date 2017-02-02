/// <reference path="ForceField.ts" />

namespace MolecularMechanics {
    export class System {
        private dihedralsCos: ForceField.DihedralCosine[];
        private dihedralsRB: ForceField.DihedralRB[];
        private angles: ForceField.Angle[];
        private bonds: ForceField.Bond[];
        private atoms: ForceField.Atom[];
 
        constructor() {
            this.dihedralsCos = [];
            this.dihedralsRB = [];
            this.angles = [];
            this.bonds = [];
            this.atoms = [];
        }

        addBond(a1: number, a2: number, kf: number, r0: number): ForceField.Bond {
            var item = new ForceField.Bond(a1, a2, kf, r0);
            this.bonds.push(item);
            return item;
        }

        addAngle(a1: number, a2: number, a3: number, kf: number, r0: number): ForceField.Angle {
        var item = new ForceField.Angle(a1, a2, a3, kf, r0);
            this.angles.push(item);
            return item;
        }

        addCosineDihedral(a1: number, a2: number, a3: number, a4: number, 
                kf: number, pn: number, ph: number): ForceField.DihedralCosine {
            var item = new ForceField.DihedralCosine(a1, a2, a3, a4, kf, pn, ph);
            this.dihedralsCos.push(item);
            return item;
        }

        addRBDihedral(a1: number, a2: number, a3: number, a4: number,
                        c0: number, c1: number, c2: number, c3: number, c4: number, c5: number): ForceField.DihedralRB {
            var item = new ForceField.DihedralRB(a1, a2, a3, a4, c0, c1, c2, c3, c4, c5);
            this.dihedralsRB.push(item);
            return item;
        }
        
        addAtom(id: number, name: string, sigma: number, epsilon: number, charge: number): ForceField.Atom {
            var item = new ForceField.Atom(id, name, sigma, epsilon, charge);
            this.atoms.push(item);
            return item;
        }

        eval(coords: Vector[], forces?: Vector[]): number {
            var Ebonded = 0;
            var Enonbonded;
            Ebonded += ForceField.Bond.eval(this.bonds, coords, forces);
            Ebonded += ForceField.Angle.eval(this.angles, coords, forces);
            Ebonded += ForceField.DihedralRB.eval(this.dihedralsRB, coords, forces);
            Ebonded += ForceField.DihedralCosine.eval(this.dihedralsCos, coords, forces);
            Enonbonded = ForceField.Atom.eval(this.atoms, coords, forces);
            return Ebonded + Enonbonded.coulomb + Enonbonded.vdw;
        }
    
    }
}
