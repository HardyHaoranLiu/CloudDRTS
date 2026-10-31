# Third-party notices

ANDES is developed by Hantao Cui and contributors, CURENT. Upstream source: https://github.com/CURENT/andes . The supplied ANDES 1.7.5 source distribution and its case data are redistributed under the GNU General Public License version 3; the unmodified license is included as `external_run/historical/LICENSE_ANDES_GPL3.txt` and inside the source archive. ANDES is provided without warranty as detailed in that license.

The official input `external_run/case/kundur_full.xlsx` is byte-identical to the file inside the bundled ANDES 1.7.5 source distribution. Its SHA-256 is `f725e03ba12d8207616f68acdd606bbd35e7c4a68f13e66d7db43925adac2ed8`.

The source archive `external_run/historical/andes-1.7.5.tar.gz` was obtained from the official PyPI distribution, verified against the PyPI SHA-256 `878e187db89484f6e778ea19f31a8830440f7bef28d4a59e115847d54092b678`. `SOURCE.json` records its release commit. No vendor equations were patched. The serial CPU scheduling shim and the explicit event initialization adapter are separate launch/application code and are disclosed in the report.

ANDES 1.9.3 is installed from PyPI and is not bundled as an installed environment. Its upstream source is available at https://pypi.org/project/andes/1.9.3/ and https://github.com/CURENT/andes/tree/v1.9.3 . The exact dependency versions used in this execution are recorded in `external_run/requirements_lock.txt`; each dependency retains its upstream license.

Only a small, attributed numerical subset from the official ANDES documentation is transcribed for verification. The full documentation is not redistributed. CoLib is cited for the parameter and load conventions, without redistribution of its page. STEPSS results are mentioned as supplementary published references; STEPSS software, result files, and data are not redistributed or claimed to have been executed.
