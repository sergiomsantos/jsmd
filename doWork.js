importScripts('ff.js');
self.addEventListener('message', function(e) {
  var data = e.data;
  switch (data.cmd) {
    case 'start':
      self.postMessage('WORKER STARTED: ' + data.msg);
      // molecule = data
      integrator.run(mm);
      break;
    case 'stop':
      self.postMessage('WORKER STOPPED');
      //self.close(); // Terminates the worker.
      // console.log(MolecularMechanics);
      break;
    default:
      self.postMessage('Unknown command: ' + data.msg);
  };
}, false);

var molecule;

var ffield = {
  atoms: 
   [ { charge: -0.507192,
       epsilon: 0.71128,
       id: 1,
       name: 'N1',
       sigma: 3.25 },
     { charge: -0.507192,
       epsilon: 0.71128,
       id: 2,
       name: 'N2',
       sigma: 3.25 },
     { charge: -0.507192,
       epsilon: 0.71128,
       id: 3,
       name: 'N3',
       sigma: 3.25 } ],
  bonds: 
   [ { a1: 0, a2: 1, kf: 3671.9, r0: 1.371 },
     { a1: 1, a2: 2, kf: 3671.9, r0: 1.371 },
     { a1: 2, a2: 0, kf: 3671.9, r0: 1.371 } ] }

var mm = new MolecularMechanics.System();

mm.loadParamsFromJSON(ffield);

mm.loadState({
  coords: [
    [0.0, 0.0, 0.0],
    [1.5, 0.0, 0.0],
    [1.5, 1.5, 0.0]
    ]
});

var fcn = function(sys, step) {
  // var crd = sys.getCoordinates();
  // self.postMessage({crd: crd, step: step});
  self.postMessage(step);
  
}

var integrator = new Driver.LeapfrogIntegrator({dt: 0.001, nsteps: 1000000});
//var thermostat = new Thermostat.Berendsen({temperature:300, tcoupl: 0.5, dt:0.0001});
var thermostat = new Thermostat.VRescale({temperature:300});
integrator.setCallback(fcn,500);
integrator.setThermostat(thermostat);
thermostat.assignVelocities(mm);
// var Ti = thermostat.getInstantTemperature(mm.getVelocities(), mm.getAtoms().map(m=>m.mass), mm.countDegreesOfFreedom());
