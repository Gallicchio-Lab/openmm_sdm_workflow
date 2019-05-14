import sys
import re
import fileinput

#usage: fpt_bound_unbound.py r*/repl.cycle.totE.potE.temp.lambda.ebind.lambda1.lambda2.alpha.u0.w0.dat

t_binding_events = 0
t_unbinding_events = 0
nfiles = len(sys.argv[1:])

for fid in range(0,nfiles):
    filename = sys.argv[fid+1]
    f = open(filename,"r")
    string = ""
    for line in f:
        line.rsplit()
        items = line.split()
        time = int(items[1])
        binding_energy = float(items[6])
        if binding_energy < -10.0:
            string += "b"
        elif binding_energy > 30.0:
            string += "u"
        else:
            string += "n"
    f.close()
    
    matches = re.findall("u[^ub]+b",string)
    binding_events = len(matches)
    fpts = [len(k) for k in matches]
    
    matches = re.findall("b[^ub]+u",string)
    unbinding_events = len(matches)
    fpts = [len(k) for k in matches]
    
    print("%s binding_events: %s unbinding_events: %d" % (filename, binding_events, unbinding_events))
    
    t_binding_events += binding_events
    t_unbinding_events += unbinding_events
    
print("Totals:  %d %d" % (t_binding_events, t_unbinding_events))
#if nevents > 0:
#    print(fpts)
