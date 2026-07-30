import logging

from byted_pysisyphus.dynamics.colvars import get_colvar
from byted_pysisyphus.dynamics.Gaussian import Gaussian
from byted_pysisyphus.dynamics.helpers import get_mb_velocities_for_geom
from byted_pysisyphus.dynamics.mdp import mdp
from byted_pysisyphus.dynamics.rattle import rattle_closure
from byted_pysisyphus.dynamics.driver import md
from byted_pysisyphus.dynamics.wigner import get_wigner_sampler


logger = logging.getLogger("dynamics")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("dynamics.log", mode="w", delay=True)
# fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
fmt_str = "%(asctime)s - %(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
