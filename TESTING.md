# Testing and numerical reproducibility

The public research project includes a read-only GitHub Actions workflow:
[Scientific verification](https://github.com/IskanderAbdullin/Dark-Matter/actions/workflows/scientific-ci.yml).

It runs automatically on pushes to main and on pull requests, and can be triggered manually. The 2026 verified major action releases are actions/checkout@v7.0.1 and actions/setup-python@v7.1.0. Python 3.11 and 3.13 are tested. Workflow permissions are limited to contents: read, checkout does not persist authentication, and there is no use of secrets or pull_request_target.

Eight independently designed checks include: Python syntax of the full copied project; positivity and monotonicity of Bose energy-cutoff remainders; outer density exponents; two self-consistent subcritical Poisson profiles; a finite-radius critical lowered-Bose solution; Uehling-Uhlenbeck detailed-balance witness; near-critical Bose susceptibility; and an exact-vs-integrated open-boundary toy model. The original synthetic SPARC unittest suite runs separately; it does **not** download or fit the real observational catalog.

To reproduce the small suite locally, install numpy, scipy, mpmath and pytest, then run \`python -m pytest -q tests\`. To run the source's synthetic SPARC tests, use \`PYTHONPATH=code python -m unittest discover -s code -p 'test_*.py' -v\` on Linux/macOS. On Windows PowerShell use \`$env:PYTHONPATH="code"; python -m unittest discover -s code -p "test_*.py" -v\`.

Further numerical programs in continuation-02 through continuation-05 are still research prototypes; individual JSON outputs are not assumed to be observations and have differing levels of independent verification. The assertions do not prove a global nonlinear Vlasov theorem, a physically generated Bose-star distribution, cosmological transfer or universal 86-eV dark-matter mass. Those require additional physical modeling and external validation.

The migration was scoped to the dedicated science folder of the original Info repository. A pattern-based credential/contact scan did not find API keys, private keys, passwords, local home directories or email-address patterns in the 141 copied UTF-8 files. Automated pattern scans are not an absolute guarantee of absence of sensitive text. Source files were SHA-verified after migration.
