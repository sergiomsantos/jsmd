importScripts('ff.js');

var doContinue = true;

self.addEventListener('message', function(e) {
  switch (e.data.cmd) {
    case 'start':
      doContinue = true;
      integrator.run(mm);
      break;
    case 'continue':
      if (doContinue)
        integrator.run(mm);
      break;
    case 'stop':
      doContinue = false;
      break;
    default:
      break;
      // self.postMessage('Unknown command: ' + data.msg);
  };
}, false);
//self.close(); // Terminates the worker.


var mm = new MolecularMechanics.System();

var s ='{"atoms":[{"sigma":3.25,"epsilon":0.71128,"charge":-0.507192,"name":"N1","id":1},{"sigma":3.39967,"epsilon":0.35982,"charge":0.099858,"name":"C1","id":2},{"sigma":3.39967,"epsilon":0.35982,"charge":0.425456,"name":"C2","id":3},{"sigma":3.25,"epsilon":0.71128,"charge":-0.59647,"name":"N2","id":4},{"sigma":3.39967,"epsilon":0.35982,"charge":0.208498,"name":"C3","id":5},{"sigma":1.06908,"epsilon":0.06569,"charge":0.391754,"name":"H1","id":6},{"sigma":3.39967,"epsilon":0.45773,"charge":-0.597793,"name":"C4","id":7},{"sigma":2.64953,"epsilon":0.06569,"charge":0.160539,"name":"H2","id":8},{"sigma":2.64953,"epsilon":0.06569,"charge":0.160539,"name":"H3","id":9},{"sigma":2.64953,"epsilon":0.06569,"charge":0.160539,"name":"H4","id":10},{"sigma":3.39967,"epsilon":0.45773,"charge":-0.420116,"name":"C5","id":11},{"sigma":2.64953,"epsilon":0.06569,"charge":0.127425,"name":"H5","id":12},{"sigma":2.64953,"epsilon":0.06569,"charge":0.127425,"name":"H6","id":13},{"sigma":2.64953,"epsilon":0.06569,"charge":0.127425,"name":"H7","id":14},{"sigma":2.42146,"epsilon":0.06276,"charge":0.132113,"name":"H8","id":15}],"bonds":[{"a1":0,"a2":1,"kf":3671.9,"r0":1.371},{"a1":0,"a2":4,"kf":3671.9,"r0":1.371},{"a1":0,"a2":5,"kf":3402.4,"r0":1.011},{"a1":1,"a2":2,"kf":4217.5,"r0":1.371},{"a1":1,"a2":10,"kf":2822.5,"r0":1.499},{"a1":2,"a2":3,"kf":3611.6,"r0":1.376},{"a1":2,"a2":6,"kf":2822.5,"r0":1.499},{"a1":3,"a2":4,"kf":4138.8,"r0":1.335},{"a1":4,"a2":14,"kf":2979,"r0":1.079},{"a1":6,"a2":7,"kf":2822.5,"r0":1.092},{"a1":6,"a2":8,"kf":2822.5,"r0":1.092},{"a1":6,"a2":9,"kf":2822.5,"r0":1.092},{"a1":10,"a2":11,"kf":2822.5,"r0":1.092},{"a1":10,"a2":12,"kf":2822.5,"r0":1.092},{"a1":10,"a2":13,"kf":2822.5,"r0":1.092}],"angles":[{"a1":0,"a2":1,"a3":2,"kf":610.03,"r0":109.42},{"a1":0,"a2":1,"a3":10,"kf":548.1,"r0":122.8},{"a1":0,"a2":4,"a3":3,"kf":625.93,"r0":112.02},{"a1":0,"a2":4,"a3":14,"kf":416.73,"r0":122.1},{"a1":1,"a2":0,"a3":4,"kf":576.56,"r0":109.9},{"a1":1,"a2":0,"a3":5,"kf":394.97,"r0":124.66},{"a1":1,"a2":2,"a3":3,"kf":594.13,"r0":114.98},{"a1":1,"a2":2,"a3":6,"kf":542.25,"r0":119.45},{"a1":1,"a2":10,"a3":11,"kf":394.97,"r0":110.86},{"a1":1,"a2":10,"a3":12,"kf":394.97,"r0":110.86},{"a1":1,"a2":10,"a3":13,"kf":394.97,"r0":110.86},{"a1":2,"a2":1,"a3":10,"kf":542.25,"r0":119.45},{"a1":2,"a2":3,"a3":4,"kf":589.94,"r0":107.47},{"a1":2,"a2":6,"a3":7,"kf":394.97,"r0":110.86},{"a1":2,"a2":6,"a3":8,"kf":394.97,"r0":110.86},{"a1":2,"a2":6,"a3":9,"kf":394.97,"r0":110.86},{"a1":3,"a2":2,"a3":6,"kf":550.61,"r0":121.44},{"a1":3,"a2":4,"a3":14,"kf":419.24,"r0":125.38},{"a1":4,"a2":0,"a3":5,"kf":394.97,"r0":124.66},{"a1":7,"a2":6,"a3":8,"kf":329.7,"r0":108.35},{"a1":7,"a2":6,"a3":9,"kf":329.7,"r0":108.35},{"a1":8,"a2":6,"a3":9,"kf":329.7,"r0":108.35},{"a1":11,"a2":10,"a3":12,"kf":329.7,"r0":108.35},{"a1":11,"a2":10,"a3":13,"kf":329.7,"r0":108.35},{"a1":12,"a2":10,"a3":13,"kf":329.7,"r0":108.35}],"dihedralsRB":[{"a1":0,"a2":1,"a3":2,"a4":3,"c0":33.472,"c1":0,"c2":-33.472,"c3":0,"c4":0,"c5":0},{"a1":0,"a2":1,"a3":2,"a4":6,"c0":33.472,"c1":0,"c2":-33.472,"c3":0,"c4":0,"c5":0},{"a1":0,"a2":1,"a3":10,"a4":11,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":0,"a2":1,"a3":10,"a4":12,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":0,"a2":1,"a3":10,"a4":13,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":0,"a2":4,"a3":3,"a4":2,"c0":39.748,"c1":0,"c2":-39.748,"c3":0,"c4":0,"c5":0},{"a1":1,"a2":0,"a3":4,"a4":3,"c0":14.2256,"c1":0,"c2":-14.2256,"c3":0,"c4":0,"c5":0},{"a1":1,"a2":0,"a3":4,"a4":14,"c0":14.2256,"c1":0,"c2":-14.2256,"c3":0,"c4":0,"c5":0},{"a1":1,"a2":2,"a3":3,"a4":4,"c0":39.748,"c1":0,"c2":-39.748,"c3":0,"c4":0,"c5":0},{"a1":1,"a2":2,"a3":6,"a4":7,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":1,"a2":2,"a3":6,"a4":8,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":1,"a2":2,"a3":6,"a4":9,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":2,"a2":1,"a3":10,"a4":11,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":2,"a2":1,"a3":10,"a4":12,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":2,"a2":1,"a3":10,"a4":13,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":2,"a2":3,"a3":4,"a4":14,"c0":39.748,"c1":0,"c2":-39.748,"c3":0,"c4":0,"c5":0},{"a1":3,"a2":2,"a3":1,"a4":10,"c0":33.472,"c1":0,"c2":-33.472,"c3":0,"c4":0,"c5":0},{"a1":3,"a2":2,"a3":6,"a4":7,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":3,"a2":2,"a3":6,"a4":8,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":3,"a2":2,"a3":6,"a4":9,"c0":0,"c1":0,"c2":0,"c3":0,"c4":0,"c5":0},{"a1":4,"a2":0,"a3":1,"a4":2,"c0":14.2256,"c1":0,"c2":-14.2256,"c3":0,"c4":0,"c5":0},{"a1":4,"a2":0,"a3":1,"a4":10,"c0":14.2256,"c1":0,"c2":-14.2256,"c3":0,"c4":0,"c5":0},{"a1":4,"a2":3,"a3":2,"a4":6,"c0":39.748,"c1":0,"c2":-39.748,"c3":0,"c4":0,"c5":0},{"a1":5,"a2":0,"a3":1,"a4":2,"c0":14.2256,"c1":0,"c2":-14.2256,"c3":0,"c4":0,"c5":0},{"a1":5,"a2":0,"a3":1,"a4":10,"c0":14.2256,"c1":0,"c2":-14.2256,"c3":0,"c4":0,"c5":0},{"a1":5,"a2":0,"a3":4,"a4":3,"c0":14.2256,"c1":0,"c2":-14.2256,"c3":0,"c4":0,"c5":0},{"a1":5,"a2":0,"a3":4,"a4":14,"c0":14.2256,"c1":0,"c2":-14.2256,"c3":0,"c4":0,"c5":0},{"a1":6,"a2":2,"a3":1,"a4":10,"c0":33.472,"c1":0,"c2":-33.472,"c3":0,"c4":0,"c5":0}],"dihedralsCos":[{"a1":0,"a2":1,"a3":2,"a4":10,"kf":4.6024,"ph":180,"pn":2},{"a1":5,"a2":0,"a3":4,"a4":1,"kf":4.6024,"ph":180,"pn":2},{"a1":6,"a2":1,"a3":2,"a4":3,"kf":4.6024,"ph":180,"pn":2},{"a1":14,"a2":0,"a3":4,"a4":3,"kf":4.6024,"ph":180,"pn":2}]}';

ffield = JSON.parse(s);
mm.loadParamsFromJSON(ffield);

mm.loadState({
  coords: [
    [ 1.988000 , 3.489000 , 0.941000 ],
    [ 3.168000 , 4.142000 , 1.181000 ],
    [ 4.140000 , 3.183000 , 1.111000 ],
    [ 3.586000 , 1.932000 , 0.979000 ],
    [ 2.271000 , 2.147000 , 0.865000 ],
    [ 1.170000 , 3.922000 , 0.533000 ],
    [ 5.461000 , 3.403000 , 1.785000 ],
    [ 5.757000 , 4.452000 , 1.718000 ],
    [ 6.240000 , 2.799000 , 1.314000 ],
    [ 5.408000 , 3.130000 , 2.841000 ],
    [ 3.456000 , 5.528000 , 0.687000 ],
    [ 2.548000 , 5.997000 , 0.300000 ],
    [ 4.198000 , 5.506000 ,-0.114000 ],
    [ 3.844000 , 6.152000 , 1.494000 ],
    [ 1.503000 , 1.393000 , 0.951000 ]
  ]
});

var callback = function(sys, step) {
  var crd = sys.getCoordinates();
  self.postMessage({crd: crd, step: step, n:VectorPool.counter});
}

var integrator = new Driver.LeapfrogIntegrator({dt: 0.001, nsteps: 100});
//var thermostat = new Thermostat.Berendsen({temperature:300, tcoupl: 0.5, dt:0.001});
var thermostat = new Thermostat.VRescale({temperature:300});
integrator.setCallback(callback, 500);
integrator.setThermostat(thermostat);
thermostat.assignVelocities(mm);
