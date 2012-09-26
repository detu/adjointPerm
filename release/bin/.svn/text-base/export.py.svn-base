#!/usr/bin/env python

# $Id: export.py 1944 2009-03-30 13:25:22Z bska $

import os, errno

# Fine scale, two-phase incompressible flow support:
def finescale():
    src = {'.'              : ['startup.m',
                               'info.xml',
                               'AUTHORS',
                               'COPYING',
                              ],
           'doc'            : ['get_start.html',
                               'helptoc.xml',
                               'top_example.html',
                              ],
           'examples'       : ['demos.xml'],
           'examples/1ph'   : ['gravityColumn.m',
                               'realField1phExample.m',
                               'simpleBC.m',
                               'simpleCornerPointExample.m',
                               'simpleSRCandBC.m',
                               'simpleTPFA.m',
                               'simpleWellExample.m',
                              ],
           'examples/2ph'   : ['realField2phExample.m',
                               'simple2phWellExample.m',
                              ],
           'examples/grids' : ['cornerPointModelExample.m',
                               #'GSmodel.grdecl',
                               'realFieldModelExample.m',
                              ],
           'gridprocessing' : ['boundaryFaceIndices.m',
                               'boundaryFaces.m',
                               'buildCornerPtNodes.m',
                               'cartGrid.m',
                               'computeGeometry.m',
                               'grid_structure.m',
                               'processGRDECL.m',
                               'removeCells.m',
                               'tensorGrid.m',
                              ],
           'gridprocessing/testgrids' :
                              ['simpleGrdecl.m'],
           'params'         : ['gravity.m'],
           'params/fluid'   : ['initFluid.m',
                               'initSimpleFluid.m',
                               'initSingleFluid.m',
                              ],
           'params/rock'    : ['logNormLayers.m',
                               'poreVolume.m',
                              ],
           'params/wells_and_bc' :
                              ['addBC.m',
                               'addSource.m',
                               'addWell.m',
                               'fluxside.m',
                               'pside.m',
                               'verticalWell.m',
                              ],
           'plotting'       : ['plotCellData.m',
                               'plotFaces.m',
                               'plotGrid.m',
                               'plotWell.m',
                              ],
           'solvers'        : ['initResSol.m',
                               'initWellSol.m',
                               'computePressureRHS.m',
                              ],
           'solvers/mimetic': ['computeMimeticIP.m',
                               'assembleWellSystem.m',
                               'packageWellSol.m',
                               'solveIncompFlow.m',
                               'unpackWellSystemComponents.m',
                              ],
           'solvers/mimetic/utils' :
                              ['mixedSymm.m',
                               'schurComplementSymm.m',
                               'tpfSymm.m',
                              ],
           'solvers/tpfa'   :
                              ['computeTrans.m',
                               'incompTPFA.m',
                              ],
           'solvers/transport':
                              ['explicitTransport.m',
                               'findFaceMobMat.m',
                               'implicitTransport.m',
                               'initFaceMob.m',
                               'initTransport.m',
                               'twophaseUpwBEGrav.m',
                               'twophaseUpwBE.m',
                               'twophaseUpwFEGrav.m',
                               'twophaseUpwFE.m',
                              ],
           'utils'          : ['cellFlux2faceFlux.m',
                               'dispif.m',
                               'faceFlux2cellFlux.m',
                               'mcolon.m',
                               'merge_options.m',
                               'msgid.m',
                               'rldecode.m',
                               'rlencode.m',
                               'ROOTDIR.m',
                               'ticif.m',
                               'tocif.m',
                              ],
           'utils/inout'    : ['readGRDECL.m',
                               'readVector.m',
                               'writeGRDECL.m',
                               #'readWellKW.m',
                              ],
           'utils/units'    : ['atm.m',
                               'barsa.m',
                               'centi.m',
                               'convertFrom.m',
                               'convertTo.m',
                               'darcy.m',
                               'day.m',
                               'deci.m',
                               'ft.m',
                               'giga.m',
                               'hour.m',
                               'inch.m',
                               'kilo.m',
                               'lbf.m',
                               'mega.m',
                               'meter.m',
                               'micro.m',
                               'milli.m',
                               'minute.m',
                               'Pascal.m',
                               'poise.m',
                               'psia.m',
                               'second.m',
                               'stb.m',
                               'year.m',
                              ],
          }

    return src

# Multiscale, two-phase incompressible flow support:
def multiscale():
    src = {
           'examples/1ph'   : ['gravityColumnMS.m',
                               'simpleBCMS.m',
                               'simpleCornerPointExampleMS.m',
                               'simpleSRCandBCMS.m',
                               'simpleWellExampleMS.m',
                               'simpleWellOverlap.m',
                              ],
           'solvers/msmfem' : ['assignBasisFuncs.m',
                               'compressPartition.m',
                               'evalBasisFunc.m',
                               'evalBasisSource.m',
                               'evalWellBasis.m',
                               'generateCoarseGrid.m',
                               'generateCoarseSystem.m',
                               'generateCoarseWellSystem.m',
                               'partitionCartGrid.m',
                               'partitionLayers.m',
                               'partitionUI.m',
                               'processPartition.m',
                               'solveIncompFlowMS.m',
                               'subFaces.m',
                               'unpackWellSystemComponentsMS.m',
                              ],
           'solvers/msmfem/utils/' :
                              ['extractBF.m',
                               'extractWellBF.m',
                               'msMatrixStructure.m',
                              ],
           'plotting'       : ['outlineCoarseGrid.m',
                               'plotBlockAndNeighbors.m'
                              ],
          }

    return src

def mkdirs(dirs, mode = 0777):
    for d in dirs:
        try:
            os.makedirs(d, mode)
        except OSError, err:
            # Reraise error if not concerned with existing dir
            if err.errno != errno.EEXIST or not os.path.isdir(d):
                raise


def export(srcdir, destdir, tree, fs, op):
    for dir in tree.keys():
        for f in tree[dir]:
            p = dir     + fs + f
            s = srcdir  + fs + p
            d = destdir + fs + p

            try:
                op(s, d)
            except Exception, e:
                if e.errno != errno.EEXIST:
                    print 'Unexpected error: ', e, ' during export()'
                    raise


def main():
    fs     = os.path.sep
    curdir = os.getcwd()
    dstdir = curdir + fs + 'export'

    pkg    = [finescale, multiscale]
    #pkg = [finescale]
    for p in pkg:
        tree = p()

        try:
            mkdirs([dstdir + fs + x for x in tree.keys()])
            export(curdir + fs + '..', dstdir, tree, fs, os.symlink)

        except Exception, e:
            print 'Failed to export package files due to ', e

    else:
        print 'Package successfully exported (by \'%s\') to\n\n\t%s' % \
              ('ln -s', dstdir)


if __name__ == '__main__':
    main()
