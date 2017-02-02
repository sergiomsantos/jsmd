from SimpleXMLRPCServer import SimpleXMLRPCServer
from threading import Thread
from pymol import cmd, util
import cPickle
import string
import os
import re

def pymol_startup () :
  print "Loading PyMOL XML-RPC extensions"
  server = pymol_interface()

class pymol_xmlrpc_server (SimpleXMLRPCServer) :
  def __init__ (self, addr, interface) :
    self.interface = interface
    SimpleXMLRPCServer.__init__(self, addr, logRequests=0)

  def _dispatch (self, method, params) :
    if not self.interface.enable_xmlrpc :
      return -1
    result = -1
    func = None
    if hasattr(self.interface, method) :
      func = getattr(self.interface, method)
    elif hasattr(cmd, method) :
      func = getattr(cmd, method)
    if not callable(func) :
      print "%s is not a callable object" % method
    else :
      result = func(*params)
      if result is None :
        result = -1
    return result

class pymol_interface (object) :
  def __init__ (self) :
    self.enable_xmlrpc = True
    # the port can be set via an environment variable - although it could just as easily be passed to __init__
    port = string.atoi(os.environ.get("PYMOL_XMLRPC_PORT", "9123"))
    self._server = pymol_xmlrpc_server(("localhost", port), self)
    t = threading.Thread(target=self._server.serve_forever)
    t.setDaemon(1)
    t.start()
    print "Started XML-RPC server on port %d" % port

  # def close_all (self) :
  #   cmd.delete("*")

  # def disable_all (self) :
  #   cmd.disable("*")

  # def load_current_model_and_maps (self,
  #     pdb_file,
  #     fwt_map,
  #     delfwt_map) :
  #   model_name = os.path.basename(os.path.splitext(pdb_file)[0])
  #   cmd.load(pdb_file, model_name, state=1)
  #   cmd.load(fwt_map, "2fofc", state=1, format="ccp4")
  #   cmd.load(delfwt_map, "fofc", state=1, format="ccp4")
  #   cmd.isomesh("m1", "2fofc", 1.0, model_name, 5.0)
  #   cmd.isomesh("m2", "fofc", 3.0, model_name, 5.0)
  #   cmd.isomesh("m3", "fofc", -3.0, model_name, 5.0)
  #   cmd.color("marine", "m1")
  #   cmd.color("green", "m2")
  #   cmd.color("red", "m3")

  # def setup_roving (self, f_map="2fofc", diff_map="fofc") :
  #   cmd.set("roving_detail", 20)
  #   cmd.set("roving_isomesh", 20)
  #   cmd.set("roving_origin", 1)
  #   cmd.set("roving_sticks", 0)
  #   cmd.set("roving_ribbon", 0)
  #   cmd.set("roving_lines", 0)
  #   cmd.set("roving_spheres", 0)
  #   cmd.set("roving_nb_spheres", 0)
  #   cmd.set("roving_polar_contacts", 0)
  #   cmd.set("roving_polar_cutoff", 0)
  #   cmd.set("stick_radius", 0.12)
  #   cmd.set("roving_map1_name", f_map)
  #   cmd.set("roving_map1_level", 1)
  #   cmd.set("roving_map2_name", diff_map)
  #   cmd.set("roving_map3_name", diff_map)
  #   cmd.set("roving_map2_level", 3)
  #   cmd.set("roving_map3_level", -3)
  #   cmd.refresh()

  # def refresh_maps (self, f_map="2fofc", diff_map="fofc") :
  #   cmd.isomesh("rov_m1", f_map, 1.0, "center", 20)
  #   cmd.isomesh("rov_m2", diff_map, 3.0, "center", 20)
  #   cmd.isomesh("rov_m3", diff_map, -3.0, "center", 20)
  #   cmd.color("skyblue", "rov_m1")
  #   cmd.color("green", "rov_m2")
  #   cmd.color("red", "rov_m3")

  # def show_selection (self, selection_string, zoom=True): #"True") :
  #   cmd.hide("sticks")
  #   util.cbay()
  #   try :
  #     cmd.show("sticks", selection_string)
  #   except Exception, e :
  #     return "Invalid selection."
  #   cmd.color("white", selection_string)
  #   if zoom :
  #     cmd.zoom(selection_string, buffer=5)

  # def recenter (self, x, y, z) :
  #   view = list(cmd.get_view())
  #   view[12] = x
  #   view[13] = y
  #   view[14] = z
  #   cmd.set_view(view)

if __name__ == "pymol" : # optional, if the module will be a command line argument to pymol
  pymol_startup()