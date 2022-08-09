# Package

version       = "0.1.0"
author        = "edoardo.giacopuzzi"
description   = "Split gene consequences from SnpEff, VEP and bcftools and select your preferred one"
license       = "MIT"
srcDir        = "src"
bin           = @["csq_selector"]
skipDirs      = @["test"]

# Dependencies

requires "nim >= 1.4.8", "hts >= 0.3.21", "zip >= 0.3.1"

