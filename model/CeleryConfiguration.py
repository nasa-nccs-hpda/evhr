from core.model.CeleryConfiguration import *

inclModules.append('evhr.model.EvhrToaCelery')
app.conf.include = inclModules
app.conf.worker_concurrency = 4
app.conf.worker_prefetch_multiplier = 1
