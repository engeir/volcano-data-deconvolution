# Changelog

## [1.0.4](https://github.com/engeir/volcano-data-deconvolution/compare/volcano-data-deconvolution-v1.0.3...volcano-data-deconvolution-v1.0.4) (2024-06-03)


### Miscellaneous

* **zenodo:** try specifying affiliation string as zenodo show it ([d8bd9ae](https://github.com/engeir/volcano-data-deconvolution/commit/d8bd9aee9efd9f12ae47b80e862a48deac0b7989))

## [1.0.3](https://github.com/engeir/volcano-data-deconvolution/compare/volcano-data-deconvolution-v1.0.2...volcano-data-deconvolution-v1.0.3) (2024-06-03)


### Bug Fixes

* **zenodo:** correctly specify zenodo metadata ([26349fb](https://github.com/engeir/volcano-data-deconvolution/commit/26349fb262e38b371e3f47ed67b31f8ac7492d0e))
* **zenodo:** only use the name field ([d05f086](https://github.com/engeir/volcano-data-deconvolution/commit/d05f0863ae17d9f4ba624f3f2a25f7be2852e783))


### Tests

* **zenodo:** test that the zenodo json schema is correct ([26349fb](https://github.com/engeir/volcano-data-deconvolution/commit/26349fb262e38b371e3f47ed67b31f8ac7492d0e))

## [1.0.2](https://github.com/engeir/volcano-data-deconvolution/compare/volcano-data-deconvolution-v1.0.1...volcano-data-deconvolution-v1.0.2) (2024-06-03)


### Documentation

* **README:** link to zenodo latest DOI ([8cccacc](https://github.com/engeir/volcano-data-deconvolution/commit/8cccacca87f4e7dde9c8f5b1371ca9ce495eda8b))
* **zenodo:** better specification of metadata name and affiliation ([8346dbb](https://github.com/engeir/volcano-data-deconvolution/commit/8346dbbaec07a0dcf4c9f614a7082a9265d1bd7b))

## [1.0.1](https://github.com/engeir/volcano-data-deconvolution/compare/volcano-data-deconvolution-v1.0.0...volcano-data-deconvolution-v1.0.1) (2024-06-03)


### Documentation

* **zenodo:** add JSON file for zenodo publication ([670431b](https://github.com/engeir/volcano-data-deconvolution/commit/670431bda8a438996091d98041093ced86a9f3f0))

## [1.0.0](https://github.com/engeir/volcano-data-deconvolution/compare/volcano-data-deconvolution-v0.7.0...volcano-data-deconvolution-v1.0.0) (2024-06-03)


### ⚠ BREAKING CHANGES

* **sim name:** change all simulations to use newer naming convention

### Features

* **cesm2 load:** optionally specify dims to average out ([875aeca](https://github.com/engeir/volcano-data-deconvolution/commit/875aecae6e3fd0745363622de2142c7e66d47809))
* clean up double waveform figure ([27f7ae1](https://github.com/engeir/volcano-data-deconvolution/commit/27f7ae1080f823c6b3203cb7120a2ee220bd6f73))
* clean up reconstruction figures ([592fce7](https://github.com/engeir/volcano-data-deconvolution/commit/592fce77e0bd32f15d9af7f42086d315e66e17e1))
* clean up reconstruction figures ([d48121f](https://github.com/engeir/volcano-data-deconvolution/commit/d48121f5954e3675334e408aa7992b30ce251817))
* **load:** add option to specify length of OB16 time series ([afe55fa](https://github.com/engeir/volcano-data-deconvolution/commit/afe55fa238f09b1aabb2baed4b850cce70ffb940))
* **load:** add time series statistics comparison class ([9baa890](https://github.com/engeir/volcano-data-deconvolution/commit/9baa8900cfa971e23c851e88802a35ff5e779aaa))
* **reconstruction:** add new spectrum comparison plot ([726843d](https://github.com/engeir/volcano-data-deconvolution/commit/726843da60f4ddfc09cd20c7af597bcf9da650c3))


### Bug Fixes

* **deconv:** forgot parenthesis that made mod 2 check pass every time ([ebff64a](https://github.com/engeir/volcano-data-deconvolution/commit/ebff64a383a3b0f0098fba5bdada1bb35d86e7d0))
* **load:** raise error for even length in deconvolution ([41718c7](https://github.com/engeir/volcano-data-deconvolution/commit/41718c7ab4c738d10ef96c81a4314ce8eea12ff6))
* **mypy:** do not call `.data` attr on numpy array object ([41718c7](https://github.com/engeir/volcano-data-deconvolution/commit/41718c7ab4c738d10ef96c81a4314ce8eea12ff6))


### Miscellaneous

* **analysis:** simple plots looking at Eurasian warming ([c4d5341](https://github.com/engeir/volcano-data-deconvolution/commit/c4d5341bfc4359edc526d0a76949a08c45841944))
* **cutoff:** update with better default for plotting cut-off ([4e7c0c5](https://github.com/engeir/volcano-data-deconvolution/commit/4e7c0c584951c6f999439d1bad1b4346b41f8644))
* **deconv:** use SO2 with decay in analysis ([ced1d03](https://github.com/engeir/volcano-data-deconvolution/commit/ced1d0310fafe8b4699a31022ea95c73858f2aa6))
* **double waveform:** finalise the plotting script ([906f719](https://github.com/engeir/volcano-data-deconvolution/commit/906f7198fa7aff00f7d66eb29c55dd02701a163f))
* **load:** use the newer tt-4sep over double-overlap ([758d154](https://github.com/engeir/volcano-data-deconvolution/commit/758d154940d475d15d1fae7a026eadea650db718))
* **mypy:** comment all false negatives ([4eef56e](https://github.com/engeir/volcano-data-deconvolution/commit/4eef56ed6462c67fcb5154467820358d4608eabb))
* **numerical soln:** create plots based on params of other sims ([0fd00e7](https://github.com/engeir/volcano-data-deconvolution/commit/0fd00e72383d7f5a204b05b750d79c1d2dad2fd0))
* **plotting:** fix time axis label and more ([e838916](https://github.com/engeir/volcano-data-deconvolution/commit/e8389162ebbb199c65f9e7f44cc9c03e5b22127b))
* save images as pdf ([a00654c](https://github.com/engeir/volcano-data-deconvolution/commit/a00654c7cde2c9aaf8d8f13d0a9c664c8beacde8))
* upload code... ([6df95bf](https://github.com/engeir/volcano-data-deconvolution/commit/6df95bf1847487a59a94d5987bb3f012be060b0b))
* use latex text style in figs, not mathrm ([78df340](https://github.com/engeir/volcano-data-deconvolution/commit/78df3409cb3cc72b72d1431a82120f5e16be4dd8))


### Styles

* fix plot styles and update sim name ([28962ce](https://github.com/engeir/volcano-data-deconvolution/commit/28962ceb9df4e42011502c95abe33beaae5b30ca))


### Code Refactoring

* **load:** slice on the properties to continue lazy loading ([4f45a37](https://github.com/engeir/volcano-data-deconvolution/commit/4f45a37e9a9a70836e1b289bb886977054901e90))
* **load:** use the deconv method of the Deconvolve object, not fppanalysis impl ([41718c7](https://github.com/engeir/volcano-data-deconvolution/commit/41718c7ab4c738d10ef96c81a4314ce8eea12ff6))
* move plot of single deconvolve object into its own private method ([b607ecb](https://github.com/engeir/volcano-data-deconvolution/commit/b607ecbd253aeec9c22ac133a081363813a6f3d7))
* remove deprecated scripts and move some nice-to-have to own dir ([308e719](https://github.com/engeir/volcano-data-deconvolution/commit/308e71968b1714c94a9133e23243a1384ca23538))
* remove old code and clean up ([2e90dc2](https://github.com/engeir/volcano-data-deconvolution/commit/2e90dc2874c645e9fe4fec3bd32175d62dd37f93))
* **sim name:** change all simulations to use newer naming convention ([930fd70](https://github.com/engeir/volcano-data-deconvolution/commit/930fd705e6587e0310e1391a9d38135aa6cca130))
* squash figures by using one x-axis ([c44ba64](https://github.com/engeir/volcano-data-deconvolution/commit/c44ba64494894ff7eaf208d4cb558324b120d284))


### Continuous Integration

* **github:** update organization name of release-please workflow ([6c1c94c](https://github.com/engeir/volcano-data-deconvolution/commit/6c1c94cce9c2e0ee9ffae45ae915ea36ebb6cac0))


### Build System

* **deps-dev:** update development dependencies ([99116ac](https://github.com/engeir/volcano-data-deconvolution/commit/99116ac50b6179c6435a28e13b9db883321ab45b))


### Documentation

* **README:** add description for entry points ([8781130](https://github.com/engeir/volcano-data-deconvolution/commit/8781130daeedfb5b289ded5b218e09cde15b113c))

## [0.7.0](https://github.com/engeir/volcano-data-deconvolution/compare/volcano-data-deconvolution-v0.6.0...volcano-data-deconvolution-v0.7.0) (2024-04-11)


### Features

* **numerical soln:** set params in json files and iterate curvefit until a fit is found ([6adffd2](https://github.com/engeir/volcano-data-deconvolution/commit/6adffd28c3cf797708f09d71e73b7ab12b684761))


### Bug Fixes

* **numeric soln:** reduce number of scaling params ([0c8f40a](https://github.com/engeir/volcano-data-deconvolution/commit/0c8f40ab304efdc8e58676c19cf49ee348fffcfa))


### Miscellaneous

* **numerical soln:** add estimate of R based on two exponential models only ([39fafe6](https://github.com/engeir/volcano-data-deconvolution/commit/39fafe6707709c5a7a8013a6a7cca91c48754129))
* **numerical soln:** save estimated parameters in dicts for later use ([9779027](https://github.com/engeir/volcano-data-deconvolution/commit/97790275707b36a36e3f4e93e8818a961de631ca))


### Styles

* **ruff:** add preview rules to lint/check and format ([4c618c6](https://github.com/engeir/volcano-data-deconvolution/commit/4c618c6d66fbee6b746dd22f9c69676be73e76bf))
* **ruff:** format using the preview settings and add code-format ([7100d09](https://github.com/engeir/volcano-data-deconvolution/commit/7100d09ac6b85f049cac6948ed8afb308dfd832b))


### Documentation

* **DOC:** improve on the linearity discussion ([39fafe6](https://github.com/engeir/volcano-data-deconvolution/commit/39fafe6707709c5a7a8013a6a7cca91c48754129))

## [0.6.0](https://github.com/engeir/volcano-data-deconvolution/compare/volcano-data-deconvolution-v0.5.0...volcano-data-deconvolution-v0.6.0) (2024-04-05)


### Features

* **analytic:** implement numerical soln to anaalytic expressions ([1ae1f78](https://github.com/engeir/volcano-data-deconvolution/commit/1ae1f781d92959eb8e8a9fbf9fa9b489d993b446))


### Bug Fixes

* **attrs:** disappearing in Deconvolution class during padding operation ([05bd23a](https://github.com/engeir/volcano-data-deconvolution/commit/05bd23af7761ee60424fb2ff0d3a5a728da28a20))


### Miscellaneous

* **analytic soln:** start implementation of numerical soln to analytic expression ([97c06f3](https://github.com/engeir/volcano-data-deconvolution/commit/97c06f30a33185e346306f8190250aef48f5b1be))

## [0.5.0](https://github.com/engeir/volcano-data-deconvolution/compare/volcano-data-deconvolution-v0.4.0...volcano-data-deconvolution-v0.5.0) (2024-04-04)


### Features

* **cut off:** implement class to cut response functions and recreate with noise ([cd4cef6](https://github.com/engeir/volcano-data-deconvolution/commit/cd4cef6fbb8c1abeab8ea6208d7fd132d7ef29f7))
* **cut off:** implement cut off analysis ([3afc16d](https://github.com/engeir/volcano-data-deconvolution/commit/3afc16d32e7959b9d630c6a4ef32b3562e0a3ddf))
* **cut off:** implement cut off test of response functions ([3775fff](https://github.com/engeir/volcano-data-deconvolution/commit/3775fff83a5810af4617dc4c81c49b293c832442))
* **cut off:** scale single response to fit with any(?) super position ([9761f46](https://github.com/engeir/volcano-data-deconvolution/commit/9761f46c59558c72743b82691c1d4ef93ed2a13c))
* double waveform, T2RF response plots, and more ([715e608](https://github.com/engeir/volcano-data-deconvolution/commit/715e6087fe293ea0f905545edf5ff3b265d4bc3e))
* **load:** add RF control data property to data and deconvolve classes ([696f789](https://github.com/engeir/volcano-data-deconvolution/commit/696f7897319916950fe01413598e191542003999))
* **load:** deconvolution classes accepts custom deconvolver ([6a3b94d](https://github.com/engeir/volcano-data-deconvolution/commit/6a3b94d6a49253a72557f05ad27a855f42d5554f))
* **ob16:** create reconstruction class ([bb8f887](https://github.com/engeir/volcano-data-deconvolution/commit/bb8f8871789f953a53fd6b10a7c8de58dfeafed9))
* plots of reconstructing temp and more ([405245f](https://github.com/engeir/volcano-data-deconvolution/commit/405245f00b6b3b0d5f229b3a2011195e06c196fc))
* **tmso2:** include TMSO2 as a CESM2 output variable ([34bdaf1](https://github.com/engeir/volcano-data-deconvolution/commit/34bdaf1ee137e96b49d760033edcf7ff7801fee2))
* working first impl of parametrisation analysis ([44f1ac1](https://github.com/engeir/volcano-data-deconvolution/commit/44f1ac15f99bffef3bf72029a4178ce9869030e1))


### Bug Fixes

* **build:** remove duplicate entries in pyproject.toml ([91cb558](https://github.com/engeir/volcano-data-deconvolution/commit/91cb55870e9796cd8eac5f9c38b02fcd29c3c460))
* **README:** incorrect git clone link ([a1451ba](https://github.com/engeir/volcano-data-deconvolution/commit/a1451ba1e1d57d7655df30842df79a9fcecb6f47))


### Miscellaneous

* add tests of deconv CESM T with RF and AOD ([65ccd5e](https://github.com/engeir/volcano-data-deconvolution/commit/65ccd5ed8da1ce13c29faa7d9de7a330e3ef6353))
* **deps:** add some standard dependencies ([d283699](https://github.com/engeir/volcano-data-deconvolution/commit/d2836998284f179a4cac103f54892b6d8ac2b4b1))
* **double waveform:** reconstruct based on analytic curve fit ([34bdaf1](https://github.com/engeir/volcano-data-deconvolution/commit/34bdaf1ee137e96b49d760033edcf7ff7801fee2))
* fix similar deconvolution of OB16 for daily and monthly data ([c7dccee](https://github.com/engeir/volcano-data-deconvolution/commit/c7dcceef1a2d82b422ded3eb6356689f1cb14943))
* initial setup of the project ([57d0f1e](https://github.com/engeir/volcano-data-deconvolution/commit/57d0f1e1651a30e3407da3be11846d366966c850))
* lots more updates to plots and other unimportant things ([6a3b94d](https://github.com/engeir/volcano-data-deconvolution/commit/6a3b94d6a49253a72557f05ad27a855f42d5554f))
* **main:** release 0.1.0 ([d0d6888](https://github.com/engeir/volcano-data-deconvolution/commit/d0d6888bd50b11739551e0c327b223f5ad65b4c4))
* **main:** release 0.1.1 ([4455ea2](https://github.com/engeir/volcano-data-deconvolution/commit/4455ea2bcd63f4bb512a146de38a85391d7bdc56))
* **main:** release 0.1.2 ([#3](https://github.com/engeir/volcano-data-deconvolution/issues/3)) ([bed95b1](https://github.com/engeir/volcano-data-deconvolution/commit/bed95b16a98dfa574aa8be04b657f7db2a68b281))
* **main:** release 0.2.0 ([#4](https://github.com/engeir/volcano-data-deconvolution/issues/4)) ([866523f](https://github.com/engeir/volcano-data-deconvolution/commit/866523f329908543ba36489713a748e85b511474))
* **main:** release 0.3.0 ([#5](https://github.com/engeir/volcano-data-deconvolution/issues/5)) ([570d767](https://github.com/engeir/volcano-data-deconvolution/commit/570d76706462ad80879ce9ecc24cc30dadec0578))
* **main:** release 0.3.1 ([#6](https://github.com/engeir/volcano-data-deconvolution/issues/6)) ([0babfc3](https://github.com/engeir/volcano-data-deconvolution/commit/0babfc3d9ae5457414bce3295c2c64401ca16f9c))
* **main:** release 0.4.0 ([#7](https://github.com/engeir/volcano-data-deconvolution/issues/7)) ([c1e9cff](https://github.com/engeir/volcano-data-deconvolution/commit/c1e9cffbde2fc07d2deb81fcf80eb3483de52df8))
* push recent WIP ([dff53d4](https://github.com/engeir/volcano-data-deconvolution/commit/dff53d47b73d5aa23bd98f90893d10dc6d67e065))
* **response:** try simple functional scalings to compare simulations ([56fa769](https://github.com/engeir/volcano-data-deconvolution/commit/56fa769f39d0e88b68145a604858d384e3ecf3b2))
* small changes ([33048d8](https://github.com/engeir/volcano-data-deconvolution/commit/33048d8d67cf6c0093d541c1485bd583fa11f06c))
* **waveform:** add log scaling of SO2 ([557b42a](https://github.com/engeir/volcano-data-deconvolution/commit/557b42a9517f66de038a41d2a3b7c9bdd9cf5936))
* wip ([27bb56f](https://github.com/engeir/volcano-data-deconvolution/commit/27bb56f23c3ec179d474f457b489c94d47eb508a))


### Styles

* **format:** run updated ruff style formatting via rye ([2d933fe](https://github.com/engeir/volcano-data-deconvolution/commit/2d933fef33a1d2ab5435f51c99a6586d1d03599b))
* **pre-commit:** format TOML and YAML in pre-commit hook ([92242de](https://github.com/engeir/volcano-data-deconvolution/commit/92242de9d1f23bbe451acec317fa502ddfa80072))


### Code Refactoring

* **load:** even lazier loading of arrays in deconvolution classes ([c0b37c4](https://github.com/engeir/volcano-data-deconvolution/commit/c0b37c4bfd5c1bc6b981a79f11638faf78980a52))
* **savefig:** move figures to generated_files dir ([3b8b62d](https://github.com/engeir/volcano-data-deconvolution/commit/3b8b62da2d2bddae19e8c2d5366f995ee561bdd2))


### Continuous Integration

* **release-please:** bump to v4 and fix config issue ([bfb47e2](https://github.com/engeir/volcano-data-deconvolution/commit/bfb47e2057b6e4ffa02b19cbd9d68dbd5ac8cc94))


### Build System

* change from poetry to rye ([54afd72](https://github.com/engeir/volcano-data-deconvolution/commit/54afd72d939a4705070069f683cb00398a108c6e))


### Documentation

* add install instructions ([662b4f8](https://github.com/engeir/volcano-data-deconvolution/commit/662b4f899d83d1699501906c018d0b89fedaed40))
* change descriptions to use rye rather than poetry ([d53dd43](https://github.com/engeir/volcano-data-deconvolution/commit/d53dd43c87bd76abf4eaa1789c0d6882ac579cb2))
* **fix:** !INFO not available as highlight, use !NOTE ([e9ee60f](https://github.com/engeir/volcano-data-deconvolution/commit/e9ee60f2e80dd2b2b9991cadac188b7ea34a0a54))
* start taking some notes for paper ([6a3b94d](https://github.com/engeir/volcano-data-deconvolution/commit/6a3b94d6a49253a72557f05ad27a855f42d5554f))
* update DOC ([c537a50](https://github.com/engeir/volcano-data-deconvolution/commit/c537a507c6cc76967cce18f7dfcb03fefbbaf9f2))

## [0.4.0](https://github.com/engeir/volcano-data-deconvolution/compare/v0.3.1...v0.4.0) (2024-03-26)


### Features

* **cut off:** implement cut off analysis ([3afc16d](https://github.com/engeir/volcano-data-deconvolution/commit/3afc16d32e7959b9d630c6a4ef32b3562e0a3ddf))
* **cut off:** scale single response to fit with any(?) super position ([9761f46](https://github.com/engeir/volcano-data-deconvolution/commit/9761f46c59558c72743b82691c1d4ef93ed2a13c))

## [0.3.1](https://github.com/engeir/volcano-data-deconvolution/compare/v0.3.0...v0.3.1) (2024-03-22)


### Bug Fixes

* **build:** remove duplicate entries in pyproject.toml ([91cb558](https://github.com/engeir/volcano-data-deconvolution/commit/91cb55870e9796cd8eac5f9c38b02fcd29c3c460))


### Styles

* **pre-commit:** format TOML and YAML in pre-commit hook ([92242de](https://github.com/engeir/volcano-data-deconvolution/commit/92242de9d1f23bbe451acec317fa502ddfa80072))

## [0.3.0](https://github.com/engeir/volcano-data-deconvolution/compare/v0.2.0...v0.3.0) (2024-03-22)


### Features

* **cut off:** implement class to cut response functions and recreate with noise ([cd4cef6](https://github.com/engeir/volcano-data-deconvolution/commit/cd4cef6fbb8c1abeab8ea6208d7fd132d7ef29f7))
* **load:** add RF control data property to data and deconvolve classes ([696f789](https://github.com/engeir/volcano-data-deconvolution/commit/696f7897319916950fe01413598e191542003999))


### Styles

* **format:** run updated ruff style formatting via rye ([2d933fe](https://github.com/engeir/volcano-data-deconvolution/commit/2d933fef33a1d2ab5435f51c99a6586d1d03599b))

## [0.2.0](https://github.com/engeir/volcano-data-deconvolution/compare/v0.1.2...v0.2.0) (2024-03-20)


### Features

* **cut off:** implement cut off test of response functions ([3775fff](https://github.com/engeir/volcano-data-deconvolution/commit/3775fff83a5810af4617dc4c81c49b293c832442))
* double waveform, T2RF response plots, and more ([715e608](https://github.com/engeir/volcano-data-deconvolution/commit/715e6087fe293ea0f905545edf5ff3b265d4bc3e))
* **load:** deconvolution classes accepts custom deconvolver ([6a3b94d](https://github.com/engeir/volcano-data-deconvolution/commit/6a3b94d6a49253a72557f05ad27a855f42d5554f))
* **ob16:** create reconstruction class ([bb8f887](https://github.com/engeir/volcano-data-deconvolution/commit/bb8f8871789f953a53fd6b10a7c8de58dfeafed9))
* plots of reconstructing temp and more ([405245f](https://github.com/engeir/volcano-data-deconvolution/commit/405245f00b6b3b0d5f229b3a2011195e06c196fc))
* working first impl of parametrisation analysis ([44f1ac1](https://github.com/engeir/volcano-data-deconvolution/commit/44f1ac15f99bffef3bf72029a4178ce9869030e1))


### Miscellaneous

* add tests of deconv CESM T with RF and AOD ([65ccd5e](https://github.com/engeir/volcano-data-deconvolution/commit/65ccd5ed8da1ce13c29faa7d9de7a330e3ef6353))
* **deps:** add some standard dependencies ([d283699](https://github.com/engeir/volcano-data-deconvolution/commit/d2836998284f179a4cac103f54892b6d8ac2b4b1))
* fix similar deconvolution of OB16 for daily and monthly data ([c7dccee](https://github.com/engeir/volcano-data-deconvolution/commit/c7dcceef1a2d82b422ded3eb6356689f1cb14943))
* lots more updates to plots and other unimportant things ([6a3b94d](https://github.com/engeir/volcano-data-deconvolution/commit/6a3b94d6a49253a72557f05ad27a855f42d5554f))
* push recent WIP ([dff53d4](https://github.com/engeir/volcano-data-deconvolution/commit/dff53d47b73d5aa23bd98f90893d10dc6d67e065))
* small changes ([33048d8](https://github.com/engeir/volcano-data-deconvolution/commit/33048d8d67cf6c0093d541c1485bd583fa11f06c))
* wip ([27bb56f](https://github.com/engeir/volcano-data-deconvolution/commit/27bb56f23c3ec179d474f457b489c94d47eb508a))


### Code Refactoring

* **load:** even lazier loading of arrays in deconvolution classes ([c0b37c4](https://github.com/engeir/volcano-data-deconvolution/commit/c0b37c4bfd5c1bc6b981a79f11638faf78980a52))
* **savefig:** move figures to generated_files dir ([3b8b62d](https://github.com/engeir/volcano-data-deconvolution/commit/3b8b62da2d2bddae19e8c2d5366f995ee561bdd2))


### Build System

* change from poetry to rye ([54afd72](https://github.com/engeir/volcano-data-deconvolution/commit/54afd72d939a4705070069f683cb00398a108c6e))


### Documentation

* change descriptions to use rye rather than poetry ([d53dd43](https://github.com/engeir/volcano-data-deconvolution/commit/d53dd43c87bd76abf4eaa1789c0d6882ac579cb2))
* start taking some notes for paper ([6a3b94d](https://github.com/engeir/volcano-data-deconvolution/commit/6a3b94d6a49253a72557f05ad27a855f42d5554f))
* update DOC ([c537a50](https://github.com/engeir/volcano-data-deconvolution/commit/c537a507c6cc76967cce18f7dfcb03fefbbaf9f2))

## [0.1.2](https://github.com/engeir/volcano-data-deconvolution/compare/v0.1.1...v0.1.2) (2024-01-12)


### Bug Fixes

* **README:** incorrect git clone link ([a1451ba](https://github.com/engeir/volcano-data-deconvolution/commit/a1451ba1e1d57d7655df30842df79a9fcecb6f47))

## [0.1.1](https://github.com/engeir/volcano-data-deconvolution/compare/v0.1.0...v0.1.1) (2024-01-12)


### Documentation

* add install instructions ([662b4f8](https://github.com/engeir/volcano-data-deconvolution/commit/662b4f899d83d1699501906c018d0b89fedaed40))
* **fix:** !INFO not available as highlight, use !NOTE ([e9ee60f](https://github.com/engeir/volcano-data-deconvolution/commit/e9ee60f2e80dd2b2b9991cadac188b7ea34a0a54))

## 0.1.0 (2024-01-12)


### Miscellaneous

* initial setup of the project ([57d0f1e](https://github.com/engeir/volcano-data-deconvolution/commit/57d0f1e1651a30e3407da3be11846d366966c850))
