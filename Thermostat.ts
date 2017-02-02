/// <reference path="Vectors.ts" />
/// <reference path="System.ts" />

const boltzman: number =0.0083145112119; // (kJ/(mol K))

namespace Thermostat {
  
  export interface IThermostat {
    update(velocities: Vector[], temperature: number): void;
    assignVelocities(system: MolecularMechanics.System): void;
    getInstantTemperature(velocities: Vector[], masses: number[], Ndf: number): number;
    
  }


  abstract class BaseThermostat implements IThermostat {
    protected temperature: number;
    constructor(temperature: number) {
      this.temperature = temperature;
    }

    assignVelocities(system: MolecularMechanics.System) {
      var randn = function(){
        var u = 1 - Math.random();
        var v = 1 - Math.random();
        return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0*Math.PI*v);
      }
      var masses: number[] = system.getAtoms().map(atom => atom.mass);
      var velocities = system.getVelocities();
      var f = this.temperature*boltzman;
      
      // generate random velocities withdrawn for a boltzman distribution at "temperature"
      velocities.forEach((vel,i) => vel.set(f*masses[i]*randn(), f*masses[i]*randn(), f*masses[i]*randn()));

      // instant temperature
      var Ti = this.getInstantTemperature(velocities, masses, system.countDegreesOfFreedom());
      
      // rescale velocities to yield desired initial temperature (T0)
      f = Math.sqrt(this.temperature/Ti);
      velocities.forEach(vel => vel.imult(f));

      var vcom = VectorPool.getVector();
      var vtmp = VectorPool.getVector();
      vcom.izero();
      velocities.forEach((vel, i) => {
        Vector.omult(vel, masses[i], vtmp);
        vcom.iadd(vtmp);
      });
      var totalMass = masses.reduce((a, b) => a + b, 0);
      vcom.imult(1.0/totalMass);
      velocities.forEach(v => v.isub(vcom));


    }

    getInstantTemperature(velocities: Vector[], masses: number[], Ndf: number): number {
      var Ti = masses.reduce((ekin, mass, i) => ekin+mass*velocities[i].normSq(), 0.0);
      return Ti/(Ndf*boltzman);
    }
    abstract update(velocities: Vector[], temperature: number): void; 

  }


  export class Berendsen extends BaseThermostat {
    tcoupl: number;
    dt: number;
    
    constructor(params: {temperature: number, tcoupl: number, dt: number}) {
      super(params.temperature);
      this.tcoupl = params.tcoupl;
      this.dt = params.dt;
    }

    update(velocities: Vector[], instantTemperature: number) {
      var scale = Math.sqrt(1.0+this.dt*(this.temperature/instantTemperature-1.0)/this.tcoupl);
      //scale = Math.sqrt(this.temperature/instantTemperature);
      velocities.forEach(v => v.imult(scale));
      
      // console.log(scale);
      // console.log(velocities[0]);
      
    }
  }

  export class VRescale extends BaseThermostat {
    
    constructor(params: {temperature: number}) {
      super(params.temperature);
    }

    update(velocities: Vector[], instantTemperature: number) {
      var scale = Math.sqrt(this.temperature/instantTemperature);
      velocities.forEach(v => v.imult(scale));
    }
  }

}