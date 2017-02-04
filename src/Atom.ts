/// <reference path="ForceField.ts" />

namespace ForceField {
    export class Atom {
        static MASS = {
            'H': 1,
            'C': 12,
            'N': 14,
            'O': 16,
        };
        id: number;
        name: string;
        sigma: number;
        charge: number;
        epsilon: number;
        blist: number[];
        mass: number;

        constructor(id: number, name: string, sigma: number, epsilon: number, charge: number) {
            this.sigma = sigma;
            this.epsilon = epsilon;
            this.charge = charge;
            this.id = id-1;
            this.name = name;
            this.blist = [Number.POSITIVE_INFINITY];
            this.mass = Atom.MASS[name.replace(/[0-9]+/, '')];
        }

        addBondedAtom(atomId: number) {
            if (this.blist.indexOf(atomId) < 0) {
                this.blist.push(atomId);
            }
            this.blist.sort((a,b) => {return a-b;});
        }

        static eval(atoms: Atom[], coords: Vector[], forces?: Vector[]): {coulomb: number, vdw: number} {
            var vij: Vector = VectorPool.getVector();
            var nAtoms = atoms.length;
            var bidx: number,
                bid: number;
            var atomi: Atom,
                atomj: Atom;
            var dijSq: number,
                sij: number,
                eij: number,
                fc: number,
                lj6: number,
                Ecoul: number,
                Evdw: number,
                EnergyVDW: number=0.0,
                EnergyCoul: number=0.0;
            
            for (var i=0; i<nAtoms-1; i++) {
                bidx = -1;
                atomi = atoms[i];
                // bid = atomi.blist[bidx];
                while (atomi.blist[++bidx]<i){}
                bid = atomi.blist[bidx];
                
                for (var j=i+1; j<nAtoms; j++) {

                    if (j===bid) {
                        bid = atomi.blist[++bidx];
                        continue;
                    }
                    // console.log(i,j,bid);
                    
                    Vector.osub(coords[j], coords[i], vij);
                    dijSq = vij.normSq();
                    atomj = atoms[j];
                    
                    sij = atomi.sigma + atomj.sigma;
                    eij = atomi.epsilon * atomj.epsilon;
                    lj6 = Math.pow(sij*sij/dijSq, 3);
                    
                    Evdw = 4.0 * eij * (lj6*lj6-lj6);
                    Ecoul = atomi.charge * atomj.charge / Math.sqrt(dijSq);
                    
                    if (forces) {
                        fc = (24*eij*(2*lj6*lj6-lj6) - Ecoul)/dijSq;
                        vij.imult(-fc);
                        forces[i].iadd(vij);
                        forces[j].isub(vij);
                    }
                    EnergyCoul += Ecoul;
                    EnergyVDW += Evdw;
                }
            }
            vij.free();
            // VectorPool.collect();
            return {coulomb: EnergyCoul, vdw: EnergyVDW};
        }
    }
}