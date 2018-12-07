# -*- coding: utf-8 -*-



import traceback
import signal
import time
import pdb
import sys

__version__ = "$Id: task.py 2130 2015-09-11 16:56:55Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/k2phot/task.py $"



"""
A task is a serial element of a pipeline intended to be run in parallal.
It contains a specific step of the process and helps make your pipeline
more modular. It also makes debugging error easier.

Here is an example

@task
def exampleTask(clip):
    #clip is a Clipboard object
    cache = clip['config.cachePath']
    outPath = clip['config.outputPath']

    clip['example'] = runSillyOperation(cache, out)

    #Check required key is produced
    clip['example.dummyValue']
    return clip

def runSillyOperation(cache, out):
    out = dict()
    out['dummyValue'] = 42
    return out
"""



def task(func):
    """A decorator for a task function

    The k2phot model is that you write wrapper code around your pipeline functions
    that extract values from a clipboard to pass as input arguments,
    and store the return values back in the clipboard. This wrapper code
    is called a "task".

    This decorator watches for any exceptions thrown by the task (or
    underlying code) and decides whether to fire up the debugger to
    figure out what's going on, or to fail gracefully by merely storing
    the raised exception for later debugging. This decorator considerably reduces
    the code duplication between different tasks.

    To use this code in your function, use the following import statement
    from task import task
    then simply write the text "@text" above your task. For example

    from task import task

    @task
    def exampleTask(clip):
        pass

    See an example down below.
    """

    def wrapperForTask(*args, **kwargs):
        """Decorator for a k2phot task

        Catches exceptions and either stores them or passes them to the debugger
        """
        assert(len(args) == 1)
        assert(len(kwargs) == 0)

        clip = args[0]
        if 'exception' in clip.keys():
            print( "INFO: %s not run because exception previously raised" \
                %(func.func_name))
            return clip

        if "__meta__" not in clip.keys():
            clip["__meta__"] = dict()

        debug = clip.get('config.debug', False)
        timeout_sec = clip.get('config.timeout_sec', defaultValue=0)

        #If specified, let the function timeout if it takes too long.
        #If the function takes too long, this will timeout
        if timeout_sec > 0:
            signal.signal(signal.SIGALRM, handleTimeout)
            signal.alarm(timeout_sec)


        dtype = type(clip)
        t0 = time.time()
        try:
            clip = func(*args, **kwargs)
        except SyntaxError as e:
            raise(e)
        except Exception as e:
            if debug:
                #Cancel timeout, if any
                signal.alarm(0)
                print( e)
                pdb.post_mortem(sys.exc_info()[2])
                raise e
            else:
                clip['exception'] = e
                clip['backtrace'] = traceback.format_exc()

        #Cancel any timeouts, if necessary
        signal.alarm(0)


        #Check that function returns a clipboard.
        if not isinstance(clip, dtype):
            if not isinstance(clip, dict):
                throwable = ValueError("FAIL: %s did not return a clipboard" %(func.func_name))

                if debug:
                    raise throwable
                else:
                    print( throwable)
                    clip = {'exception': throwable}
                    clip['functionName'] = func.__name__

        key = "%s-elapsedTime" %(func.__name__)
        clip["__meta__"][key] = time.time() - t0
        return clip

    wrapperForTask.__name__ = func.__name__
    return wrapperForTask



class TimeoutError(Exception):
    pass

def handleTimeout(signum, frame):
    raise TimeoutError()




    # set the timeout handler
    try:
        result = func(*args, **kwargs)
    except TimeoutError as exc:
        result = default
    finally:
        signal.alarm(0)

    return result

