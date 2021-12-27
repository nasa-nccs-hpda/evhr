from redis import exceptions
import subprocess
import sys

from core.model.CeleryConfiguration import *
from core.model.SystemCommand import SystemCommand


# -----------------------------------------------------------------------------
# class ILProcessController
#
# This class manages the lifecycle of Celery and the Redis server.
# The Python ‘with’  statement is used to manage the Celery lifecycle
# (Sec 8 PEP 343).  The 'with' statement clarifies code that previously would
# use try...finally blocks to ensure that clean-up code is executed.
# This control-flow structure has a basic structure of:
#
#   with expression [as variable]:
#       with-block
#
# The expression is evaluated and results in an object that supports the
# context management protocol (e.g.,has __enter__() and __exit__() methods).
#
# The object's __enter__() is called before with-block is executed and
# therefore can run set-up code. It also may return a value that is bound to
# the name variable, if given. Celery/Redis instantiation occurs here.
#
# After execution of the with-block is finished, the object's __exit__()
# method is called, even if the block raised an exception, and can therefore
# run clean-up code.  Celery/Redis shutdown and cleanup occurs here.  The
# application business logic is what is implemented in the with-block between
# the enter/exit calls.  The remaining structure of the IL client application
# remains unchanged.
#
# Standard Python system utilities (i.e., os.system(), subprocess.run(),
# kill(), pkill()) are used to invoke parameterized commands for:
# a) starting the Redis server
# b) registering the Celery workers, and
# c) shutting down the Redis server and Celery workers.
#
# Note that the Redis server acts as the backend that maintains state for
# Celery transactions.  This could have been implemented with brokers like
# RabbitMQ or Amazon SQS; however, Redis presented installation simplicity.
# -----------------------------------------------------------------------------
class ILProcessController():

    backendProcessId = 0

    # -------------------------------------------------------------------------
    # __enter__
    #
    # Start Redis server and Celery workers
    # -------------------------------------------------------------------------
    def __enter__(self):

        print('In ILProcessController.__enter__()')

        try:

            # Retrieve configuration port - default to 6380
            _backendPort = app.conf.get(IL_PORT)
            _backendPort = '6379' if _backendPort is None else _backendPort

            # Start the Celery Server
            ILProcessController.backendProcessId = \
                (subprocess.Popen(["/usr/local/bin/redis-server",
                                   "--protected-mode",
                                   "no",
                                   "--port",
                                   str(_backendPort)],
                                  stdout=subprocess.PIPE)).pid

            print("Redis port = ", _backendPort,
                  "ProcessId = ", ILProcessController.backendProcessId)

            # Retrieve path to configuration file
            ILProcessController.celeryConfig = 'evhr.model.' + \
                'CeleryConfiguration'

            # Retrieve concurrency level - default to max available
            _concurrency = app.conf.get(IL_CONCURRENCY)
            _concurrency = "" if _concurrency is None \
                else " --concurrency=" + _concurrency

            # Retrieve log level - default to 'info'
            _logLevel = app.conf.get(IL_LOGLEVEL)
            _logLevel = 'info' if _logLevel is None \
                else _logLevel

            # Start the Celery Workers
            _worker = "/usr/local/bin/celery -A " + \
                ILProcessController.celeryConfig + " worker " + \
                _concurrency + \
                " --loglevel=" + _logLevel + \
                " &"

            retcode = subprocess.run(_worker,
                                     shell=True,
                                     check=True,
                                     text=True)
            print(retcode)

        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)

    # -------------------------------------------------------------------------
    # __exit__
    #
    # Shutdown Redis server and Celery workers
    # -------------------------------------------------------------------------
    def __exit__(self, type, value, traceback):

        try:
            print('In ILProcessController.__exit__()',
                  ILProcessController.backendProcessId)

            # Shutdown the Celery workers
            shutdownWorkers = "/usr/bin/pkill -9 -f  " + \
                              ILProcessController.celeryConfig
            SystemCommand(shutdownWorkers, None, True)

            return True

        except exceptions.ConnectionError as inst:
            print("Connection Error ignore")
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)
        except Exception as inst:
            print(type(inst))  # the exception instance
            print(inst.args)  # arguments stored in .args
            print(inst)  # __str__ allows args to be printed directly,
