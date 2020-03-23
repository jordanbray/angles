let
  moz_overlay = import (builtins.fetchTarball https://github.com/mozilla/nixpkgs-mozilla/archive/master.tar.gz);
  nixpkgs = import <nixpkgs> { overlays = [ moz_overlay ]; };
  #rustNightlyChannel = (nixpkgs.rustChannelOf { date = "2019-01-26"; channel = "nightly"; }).rust;
  rustStableChannel = nixpkgs.latest.rustChannels.stable.rust.override {
    extensions = [
      "rust-src"
      "rls-preview"
      "clippy-preview"
      "rustfmt-preview"
    ];
  };
in
with import <nixpkgs> {}; {
  lyzsEnv = stdenv.mkDerivation {
    name = "poly";
    buildInputs = [ rustStableChannel rls pkgconfig python3 xorg.libxcb xorg.libX11 wayland wayland-protocols patchelf libGLU_combined clang llvm llvmPackages.libclang ];
  };
}

