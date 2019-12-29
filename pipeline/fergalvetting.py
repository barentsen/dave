
"""
This is a template top level script.


"""

import dave.pipeline.clipboard as clipboard

import dave.pipeline.vettingpipeline as vp
import os
import gc

def main():
    """A bare bones main program"""
    cfg = loadMyConfiguration()

    epicList = [206103150]

    for epic in epicList:
        runOne(epic, cfg)


def loadMyConfiguration():
    cfg = dict()
    cfg['debug'] = True
    cfg['campaign'] = 3
    cfg['timeout_sec'] = 120

    davePath = os.path.join(os.environ['HOME'],"daveOutput","")
    cfg['modshiftBasename'] =  davePath
    cfg['onepageBasename']  = davePath
    #Location of the place all the light curves and TPF files are stored
    cfg['dataStorePath'] = os.path.join(os.environ['HOME'],".mastio/k2")

    cfg['detrendTypes'] = ["PDC", "EVEREST", "AGP", "SFF"]

    #My front end
    tasks = """vp.serveTask vp.pickDefaultParams
                vp.makeOverviewPlotsTask""".split()
    cfg['taskList'] = tasks

    cfg['keysToIgnoreWhenSaving'] = ["serve"]
    return cfg




import multiprocessing
import contextlib
import parmap
from multiprocessing import pool
def runAll(func, iterable, config):
    """Run func over every element on iterable in parallel.

    Inputs:
    ----------
    func
	(A function) The top level function, e.g runOne(), below

    iterable
	(list, array, etc.) A list of values to operate on.

    config
	(Clipboard) A configuration clipboard.
    """

    #Turn off parallel when debugging.
    parallel = True
    if config.get('debug', False):
        parallel = False

    count = multiprocessing.cpu_count()-1

    with contextlib.closing(pool.Pool(count)) as p:
        out = parmap.map(func, iterable, config, pool=p, parallel=parallel, chunksize=5)

    return out


def runOne(k2id, config, returnClip=False):
    """Run the pipeline on a single target.

    Inputs:
    ------------
    k2id
        (int) Epic id of target to run on.

    config
        (dict) Dictionary of configuration parameters as created by, e.g
        loadMyConfiguration()

    Returns:
    ---------
    A clipboard containing the results.

    Notes:
    ---------
    Don't edit this function. The pipeline can recover gracefully from
    errors in any individual task, but an error in this function will crash
    the pipeline
    """

    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config
    clip['value'] = k2id

    #Check that all the tasks are properly defined
    print "Checking tasks exist"
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        print "Running %s" %(t)
        f = eval(t)
        clip = f(clip)

    gc.collect()
    if returnClip:
        return clip


