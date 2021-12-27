from core.model.CeleryConfiguration import *


inclModules.append('evhr.model.EvhrToaCelery')
app.conf.include = inclModules
app.conf.setdefault(IL_CONCURRENCY, '4')
