/// <reference path="ForceField.ts" />

namespace MolecularMechanics {
    
    export class System {
        private dihedralsCos: ForceField.DihedralCosine[];
        private dihedralsRB: ForceField.DihedralRB[];
        private angles: ForceField.Angle[];
        private bonds: ForceField.Bond[];
        private atoms: ForceField.Atom[];
        private x: Vector[];
        private f: Vector[];
        private v: Vector[];

        constructor() {
            this.dihedralsCos = [];
            this.dihedralsRB = [];
            this.angles = [];
            this.bonds = [];
            this.atoms = [];
        }

        setCoordinates(coords: Vector[]) {
            this.x = coords;
        }
        getCoordinates(): Vector[] {
            return this.x;
        }
        setForces(forces: Vector[]) {
            this.f = forces;
        }
        getForces(): Vector[] {
            return this.f;
        }

        setVelocities(v: Vector[]) {
            this.v = v;
        }
        getVelocities(): Vector[] {
            return this.v;
        }

        getAtoms(): ForceField.Atom[] {
            return this.atoms;
        }

        countDegreesOfFreedom(): number {
            return 3*this.atoms.length-6;
        }
        // loadState(state: {coords?: {x:number,y:number,z:number},
        //                   forces?: {x:number,y:number,z:number},
        //                   velocities?: {x:number,y:number,z:number}}) {
        loadState(state: {coords?: number[], forces?: number[], velocities?: number[]}) {
            var Natoms = this.getAtomCount();
            console.log(Natoms);
            
            this.x = Vector.getVector(Natoms);
            this.f = Vector.getVector(Natoms);
            this.v = Vector.getVector(Natoms);
            console.log(this.x);
            console.log(this.f);
            console.log(this.v);
            if (state.coords) {
                this.x.forEach((x,i) => {
                    // x.set(state.coords[i].x, state.coords[i].y, state.coords[i].z);
                    x.set(state.coords[i][0], state.coords[i][1], state.coords[i][2]);
                });
            }
            if (state.forces) {
                this.f.forEach((f,i) => {
                    // f.set(state.forces[i].x, state.forces[i].y, state.forces[i].z);
                    f.set(state.forces[i][0], state.forces[i][1], state.forces[i][2]);
                });
            }
            if (state.velocities) {
                this.v.forEach((v,i) => {
                    // v.set(state.velocities[i].x, state.velocities[i].y, state.velocities[i].z);
                    v.set(state.velocities[i][0], state.velocities[i][1], state.velocities[i][2]);
                });
            }
        }


        //loadParamsFromJSON(jsonStr: string) {
            // var obj = JSON.parse(jsonStr);
        loadParamsFromJSON(obj: {atoms: ForceField.Atom[], bonds: ForceField.Bond[],
            angles: ForceField.Angle[], dihedralsRB: ForceField.DihedralRB[],
            dihedralsCos: ForceField.DihedralCosine[]}) {
            //console.log(obj);
            
            var deg2rad = Math.PI/180.0;
            if (obj.atoms)
            obj.atoms.forEach(it => {
                console.log(this.addAtom(it.id, it.name, it.sigma, it.epsilon, it.charge));
            });

            if (obj.bonds)
            obj.bonds.forEach(it => {
                console.log(this.addBond(it.a1, it.a2, it.kf, it.r0));
            });

            if (obj.angles)
            obj.angles.forEach(it => {
                console.log(this.addAngle(it.a1, it.a2, it.a3, it.kf, deg2rad*it.r0));
            });

            if (obj.dihedralsRB)
            obj.dihedralsRB.forEach(it => {
                console.log(this.addRBDihedral(it.a1, it.a2, it.a3, it.a4,
                                               it.c0, it.c1, it.c2, it.c3, it.c4, it.c5));
            });

            if (obj.dihedralsCos)
            obj.dihedralsCos.forEach(it => {
                console.log(this.addCosineDihedral(it.a1, it.a2, it.a3, it.a4, it.kf, it.pn, deg2rad*it.ph));
            });
        }

        getAtomCount(): number {
            return this.atoms.length;
        }

        addBond(a1: number, a2: number, kf: number, r0: number): ForceField.Bond {
            var item = new ForceField.Bond(a1, a2, kf, r0);
            this.bonds.push(item);
            this.atoms[a1].addBondedAtom(a2);
            this.atoms[a2].addBondedAtom(a1);
            return item;
        }

        addAngle(a1: number, a2: number, a3: number, kf: number, r0: number): ForceField.Angle {
        var item = new ForceField.Angle(a1, a2, a3, kf, r0);
            this.angles.push(item);
            this.atoms[a1].addBondedAtom(a2);
            this.atoms[a1].addBondedAtom(a3);
            this.atoms[a2].addBondedAtom(a1);
            this.atoms[a2].addBondedAtom(a3);
            this.atoms[a3].addBondedAtom(a1);
            this.atoms[a3].addBondedAtom(a2);
            return item;
        }
        private tmp(i: number, j: number, k: number, l: number) {
            this.atoms[i].addBondedAtom(j);
            this.atoms[i].addBondedAtom(k);
            this.atoms[i].addBondedAtom(l);//<<<<<<<<
            this.atoms[j].addBondedAtom(i);
            this.atoms[j].addBondedAtom(k);
            this.atoms[j].addBondedAtom(l);
            this.atoms[k].addBondedAtom(i);
            this.atoms[k].addBondedAtom(j);
            this.atoms[k].addBondedAtom(l);
            this.atoms[l].addBondedAtom(i);//<<<<<<<<
            this.atoms[l].addBondedAtom(j);
            this.atoms[l].addBondedAtom(k);
        }
        addCosineDihedral(a1: number, a2: number, a3: number, a4: number, 
                kf: number, pn: number, ph: number): ForceField.DihedralCosine {
            var item = new ForceField.DihedralCosine(a1, a2, a3, a4, kf, pn, ph);
            this.dihedralsCos.push(item);
            this.tmp(a1,a2,a3,a4);
            return item;
        }

        addRBDihedral(a1: number, a2: number, a3: number, a4: number,
                        c0: number, c1: number, c2: number, c3: number, c4: number, c5: number): ForceField.DihedralRB {
            var item = new ForceField.DihedralRB(a1, a2, a3, a4, c0, c1, c2, c3, c4, c5);
            this.dihedralsRB.push(item);
            this.tmp(a1,a2,a3,a4);
            return item;
        }
        
        addAtom(id: number, name: string, sigma: number, epsilon: number, charge: number): ForceField.Atom {
            // fSI = 1/(4.pi.Eo) = 1389.35485 kJ mol-1 A e-2
            // sqrtfSI = sqrt(fSI);
            var sqrtfSI=37.27405062506623;
            // Lorentz-Bertelo combination rules for the 6-12 Lennard-Jones potential:
            // LJeij = (Ei*Ej)^(1/2)
            // LJsij = (Si+Sj)/2
            var item = new ForceField.Atom(id, name, 0.5*sigma, Math.sqrt(epsilon), sqrtfSI*charge);
            this.atoms.push(item);
            return item;
        }

        eval(coords: Vector[], forces?: Vector[]): number {
            var Ebonded = 0;
            var Enonbonded;
            if (forces)
                forces.forEach(force => force.izero());
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
            // return Ebonded;// + Enonbonded.coulomb + Enonbonded.vdw;
            return Ebonded + Enonbonded.coulomb + Enonbonded.vdw;
        }
    
    }
}
