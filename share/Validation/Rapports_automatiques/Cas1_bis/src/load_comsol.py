import MEDLoader as mc

def loadFile(fName):
  m = mc.MEDCouplingUMesh("mesh", 2)

  with open(fName) as f:
    l = f.readlines()
    l = [i.strip() for i in l]
    coo = []
    found = False
    for i, lin in enumerate(l):
      if lin == "# Mesh point coordinates":
        found = True
        continue
      elif not found:
        continue
      if lin == "":   break# fin des coords

      spl = lin.split()
      coo.append((float(spl[0]), float(spl[1])))
    m.setCoords(mc.DataArrayDouble(coo))

    conn = []
    found = 0
    for i, lin in enumerate(l[i:]):
      # the same for connectivity of TRI3
      if lin == "3 tri # type name":
        found = 1
        continue
      elif found == 0:
        continue
      if lin == "# Elements" and found == 1:
        found = 2
        continue
      elif found == 1:
        continue
      if lin == "":   break  # fin de la conn
      spl = lin.split()
      conn.extend([mc.NORM_TRI3, int(spl[0]), int(spl[1]), int(spl[2])])

    cI = range(0, len(conn)+1, 4)
    m.setConnectivity(mc.DataArrayInt(conn), mc.DataArrayInt(cI))
    #print dir(m)
    m.checkCoherency1()
  return m

if __name__ == "__main__":
  import sys
  m = loadFile(sys.argv[1])
  # Prepare borders
  mDesc, d, dI, r, rI = m.buildDescendingConnectivity()
  mDesc.setName("desc")
  dsi = rI.deltaShiftIndex()
  ids = dsi.getIdsEqual(1)
  print ids.getValues()
  bound = mDesc[ids]
  bound.setName("mesh")
  grp = mc.DataArrayInt(ids.getNumberOfTuples()); grp.iota(); grp.setName("bound")
  mf = mc.MEDFileUMesh()
  mf.setMeshAtLevel(0, m)
  mf.setMeshAtLevel(-1, bound)
  mf.setGroupsAtLevel(-1, grp)
  mf.write(sys.argv[2], 2)
#   mc.WriteUMesh(, m, True)
  print "OK"


