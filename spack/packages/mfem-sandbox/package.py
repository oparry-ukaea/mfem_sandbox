from spack.package import *


class MfemSandbox(CMakePackage):
    """Sandbox for testing new MFEM-based apps."""

    homepage = "https://github.com/oparry-ukaea/mfem_sandbox"
    git = "https://github.com/oparry-ukaea/mfem_sandbox"

    maintainers("oparry-ukaea")
    license("MIT", checked_by="oparry-ukaea")

    # Versions
    version("develop", branch="main")

    # Dependencies
    # N.B. static mfem doesn't link with hypre properly (hypre doesn't have a 'static' variant)
    depends_on("mfem+mpi+shared+superlu-dist cxxstd=17", type=("build", "link"))
    depends_on("mpi", type=("build", "link", "run"))

    # Variants
    # variant("variantname", default=default-val, description="desc.")

    def cmake_args(self):
        # ON/OFF definitions controlled by variants
        def_variant_map = {}
        # = {"DEF_NAME_WITHOUT_D":"variant_name"}
        variants_args = [
            self.define_from_variant(def_str, var_str)
            for def_str, var_str in def_variant_map.items()
        ]

        # Concatenate different arg types and return
        args = []
        args.extend(variants_args)

        return args
