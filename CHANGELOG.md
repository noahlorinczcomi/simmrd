# Changelog

All notable changes to this project are recorded here, organized by date.

---

## 2026-05-17

- [`df97a7b`](https://github.com/noahlorinczcomi/simmrd/commit/df97a7b) Fix `.Rbuildignore` pattern (`^cli$` → `^cli(/|$)`) so `R CMD build` correctly excludes the `cli/` subtree and its pixi environment; switch `pixi run setup` from `remotes::install_local` to `R CMD INSTALL --no-build-vignettes` to avoid a recursive source-copy that failed on broken symlinks

## 2026-05-08

- [`c6b0e9a`](https://github.com/noahlorinczcomi/simmrd/commit/c6b0e9a37e1fbb0efc5722e26d5e55e766e0a172) docs: table of outputs
- [`1a9d70a`](https://github.com/noahlorinczcomi/simmrd/commit/1a9d70a0cbfac6ccf6e23933781ee2f856c776b7) docs: cli shoutout

## 2026-04-25

- [`b10dd5c`](https://github.com/noahlorinczcomi/simmrd/commit/b10dd5cdbdd902b94bf9f2be62892d34bdca07ca) chore: code owner

## 2026-04-24

- [`8e51033`](https://github.com/noahlorinczcomi/simmrd/commit/8e51033c9c2a54a4e14242db3cc244dbffc8c610) fix: ci/cd fail
- [`81d68be`](https://github.com/noahlorinczcomi/simmrd/commit/81d68be801eb307a217aaf61771ced2ab46d0608) feat: ci/cd

## 2026-04-21

- [`397d994`](https://github.com/noahlorinczcomi/simmrd/commit/397d9946a16876aaef27c70dc9cb11d0d8a82a46) chore: refactor
- [`b578d11`](https://github.com/noahlorinczcomi/simmrd/commit/b578d11cb27394bdfb1e1657d332f91eb1f56fb2) feat: CLI

## 2026-02-16

- [`a94eac7`](https://github.com/noahlorinczcomi/simmrd/commit/a94eac7b92d6ae95a481d15ccdf46d9111debbf7) chore: update

## 2025-11-02

- [`45cbafb`](https://github.com/noahlorinczcomi/simmrd/commit/45cbafb928388f1e267390ebd075ef864b2b7c1a) Fixed bug
- [`fd72edb`](https://github.com/noahlorinczcomi/simmrd/commit/fd72edb3aced19613899f379eaa1dc78c56c15f5) Added seed option

## 2025-03-19

- [`decba75`](https://github.com/noahlorinczcomi/simmrd/commit/decba750cacb83628c9384e5ff8b38d45735a8a2) Add top-level summary paragraph to README describing what `simmrd` does and its intended use in MR method evaluation

## 2024-03-20

- [`986323e`](https://github.com/noahlorinczcomi/simmrd/commit/986323e3c913d43f1e27dc2267d9abd48fefbec4) Fix typo in RhoME description ("used my" → "used by")

## 2024-02-16

- [`2dfcf8a`](https://github.com/noahlorinczcomi/simmrd/commit/2dfcf8ac4da3d90e35f8b65a1c8734412875ab9c) Update paper citation from preprint DOI to published journal DOI (Genetic Epidemiology, gepi.22544)

## 2024-01-27

- [`3846ff9`](https://github.com/noahlorinczcomi/simmrd/commit/3846ff9f04b47401df2caefb1cb5094f351c9fcb) Add published paper citation; replace preprint-status note with citation text

## 2024-01-22

- [`cd2afe8`](https://github.com/noahlorinczcomi/simmrd/commit/cd2afe8bb746208fb7c9a8d6e50514034f7a5aa3) Rename all `simmr` references to `simmrd` throughout README (install commands, links, tutorial URLs) following package rename

## 2024-01-12

- [`f931e2c`](https://github.com/noahlorinczcomi/simmrd/commit/f931e2c1e8fdecebb715dee7bd20dc4cbed77206) Rename simmr.R to simmrd.R
- [`6b530d3`](https://github.com/noahlorinczcomi/simmrd/commit/6b530d36cf5a823d5e94778a3d2f1e6ea54cab40) Update DESCRIPTION

## 2024-01-10

- [`60ad7fb`](https://github.com/noahlorinczcomi/simmrd/commit/60ad7fbc761f4fe3c41b2aaec0ae584f8d81c53a) Update simmr.R

## 2024-01-05

- [`b99e9e8`](https://github.com/noahlorinczcomi/simmrd/commit/b99e9e86bd328920c7ef0a5b04602a5df2565502) Merge pull request [#1](https://github.com/noahlorinczcomi/simmrd/pull/1) from remlapmot/suggestions
- [`c2ebfb8`](https://github.com/noahlorinczcomi/simmrd/commit/c2ebfb8be66e38cc59cba293507f60e81f978353) Add the tutorial as a vignette
- [`339e13d`](https://github.com/noahlorinczcomi/simmrd/commit/339e13d5b6728eb2a486f40d16524c9c3b1cf9b9) Create globals.R
- [`aad32bc`](https://github.com/noahlorinczcomi/simmrd/commit/aad32bca55c80461f82fb91dd5cfaba6191640f4) Add generate_summary() example (dontrun)
- [`fb65eb9`](https://github.com/noahlorinczcomi/simmrd/commit/fb65eb9c88d26cc593058e19c09bfbb5b58b4512) Add generate_individual() example (dontrun)
- [`91f9435`](https://github.com/noahlorinczcomi/simmrd/commit/91f9435d290e40e624cf2ded720a4e3eac7399a7) Add plot_simdata() example (dontrun)
- [`042d8a5`](https://github.com/noahlorinczcomi/simmrd/commit/042d8a5f1111198dfd5bdfadd6ff9f20f233c81e) Amend description
- [`74cadab`](https://github.com/noahlorinczcomi/simmrd/commit/74cadabcd97e873c7f8ab685f3f42c0e3e279375) devtools::document()
- Helpfile and documentation fixes: argument text, `\dontrun` usage, descriptive text, heading levels, code chunk formatting, whitespace ([`9a4e446`](https://github.com/noahlorinczcomi/simmrd/commit/9a4e446afb895f7eec8069e75a0cd3157ed747a0) [`65da71f`](https://github.com/noahlorinczcomi/simmrd/commit/65da71fc25f7cb0519aaec1d1023af37eb166511) [`70f7caa`](https://github.com/noahlorinczcomi/simmrd/commit/70f7caaa3447073fd9d2d69ad68b35cace5887e4) [`33856c4`](https://github.com/noahlorinczcomi/simmrd/commit/33856c476b1897bf5994914949735bbc45c789e2) [`83a2a35`](https://github.com/noahlorinczcomi/simmrd/commit/83a2a3572e92a65b70c2d4f541c8a96ccfde07cb) [`92a5745`](https://github.com/noahlorinczcomi/simmrd/commit/92a5745df9370da01b074fc8ff4848d23c074191) [`3e5806b`](https://github.com/noahlorinczcomi/simmrd/commit/3e5806bcf95e1fd8ed2a4fa4ce11a4b516dab644) [`658554`](https://github.com/noahlorinczcomi/simmrd/commit/658554938a2b77eca8236b1302298458e4bdc733) [`bf275ee`](https://github.com/noahlorinczcomi/simmrd/commit/bf275eedad595549b0612ed892f57226d87abc40) [`f3e6c50`](https://github.com/noahlorinczcomi/simmrd/commit/f3e6c501ed59150a269405902e8b1f0e2500ff77) [`581a432`](https://github.com/noahlorinczcomi/simmrd/commit/581a432de1a7885dc525548b65c696f8633115f4) [`67d60df`](https://github.com/noahlorinczcomi/simmrd/commit/67d60df29225404ea5c746185c45f3ce95e0c00a) [`29413e3`](https://github.com/noahlorinczcomi/simmrd/commit/29413e3daf6e4495fa439c505b3c127b85cfc0d1) [`1f6fa26`](https://github.com/noahlorinczcomi/simmrd/commit/1f6fa26ed67b1e3c83c8fc18bc9cbd0dd0f1cdb3) [`88cc836`](https://github.com/noahlorinczcomi/simmrd/commit/88cc83635f0c91dbd1cb417747bc3e161ac9ebed) [`26e1423`](https://github.com/noahlorinczcomi/simmrd/commit/26e1423311642f1db861073fdf260f432e119a02) [`a27254a`](https://github.com/noahlorinczcomi/simmrd/commit/a27254a5e4df9a20a40d9362f6bd833d4b75d34f) [`0fa465d`](https://github.com/noahlorinczcomi/simmrd/commit/0fa465dc51ab32ad95887b57293adb2ceaaa339e) [`cb9dc84`](https://github.com/noahlorinczcomi/simmrd/commit/cb9dc84d5d7c62d8b44c56605a85ec0a306f7610) [`0640c66`](https://github.com/noahlorinczcomi/simmrd/commit/0640c66a5fc7296cf2b4718568c0b96f2aea1240) [`716bee5`](https://github.com/noahlorinczcomi/simmrd/commit/716bee5f76deea3cfb804610db6e2dbee9864757) [`820b59b`](https://github.com/noahlorinczcomi/simmrd/commit/820b59bd12977c971e08c0f287780fe2b93bf4eb) [`d24dd96`](https://github.com/noahlorinczcomi/simmrd/commit/d24dd9682ca8e751727b302bfe42532069ee3c25) [`4396e6c`](https://github.com/noahlorinczcomi/simmrd/commit/4396e6ce81ed02c3836bf6a3823a43d7ff933dd4) [`da52c2b`](https://github.com/noahlorinczcomi/simmrd/commit/da52c2b524a14f3afa9ae6c908f5dc164d42e58a) [`61a7119`](https://github.com/noahlorinczcomi/simmrd/commit/61a711967172f163de4633fbe3b529aa1293e358) [`f241102`](https://github.com/noahlorinczcomi/simmrd/commit/f24110257a0ae89e71375937c990f9177df7341e) [`0eb8dc8`](https://github.com/noahlorinczcomi/simmrd/commit/0eb8dc8cb1977c45b8bed63726a7de5e3d3b2dcd))

## 2024-01-04

- Add `stats::` and `graphics::` namespace qualifiers throughout ([`df3d92c`](https://github.com/noahlorinczcomi/simmrd/commit/df3d92c875e0694e568f99871fadcaaed1a23c08) [`dc8a5c7`](https://github.com/noahlorinczcomi/simmrd/commit/dc8a5c73497b56b486bead2c6c447145c2864852) [`8c04151`](https://github.com/noahlorinczcomi/simmrd/commit/8c04151c85fd3a21d04625d7d4f5186136a40d35) [`54b66ff`](https://github.com/noahlorinczcomi/simmrd/commit/54b66ff2020b08e0ee7c580aa573d3e99beef9df) [`a0b168e`](https://github.com/noahlorinczcomi/simmrd/commit/a0b168e01c7cdabadb4b0fdc767e27e98d933eb3) [`cc45122`](https://github.com/noahlorinczcomi/simmrd/commit/cc4512280d64ba84ba95bd86b41c6d1f213131d9) [`9b5cb90`](https://github.com/noahlorinczcomi/simmrd/commit/9b5cb908d65291b56ff1b31a3f0b5deea58f6634) [`281b4a6`](https://github.com/noahlorinczcomi/simmrd/commit/281b4a6b76c3fbeda6cfdad97460431d28af6d21) [`197aa4d`](https://github.com/noahlorinczcomi/simmrd/commit/197aa4d8294d8c46b782d10f2cc7a1c332fc68b0) [`9e2424f`](https://github.com/noahlorinczcomi/simmrd/commit/9e2424fc7341fabc42949d000c928cf94326d85d) [`27fc99c`](https://github.com/noahlorinczcomi/simmrd/commit/27fc99cd494f15610ad92982f48c349b21795b51) [`70bfb4c`](https://github.com/noahlorinczcomi/simmrd/commit/70bfb4c65c406e008d26939e9bbc4b17c50de31d) [`f35a106`](https://github.com/noahlorinczcomi/simmrd/commit/f35a1062d93b9bdfeed9f8708da6e906543dc841) [`c24f1be`](https://github.com/noahlorinczcomi/simmrd/commit/c24f1be22ce82354c8004dfa35fef6dedaa3fc7e) [`9f770b0`](https://github.com/noahlorinczcomi/simmrd/commit/9f770b034e30d42f11ebe537ab9a73ad6d9cd33a) [`8e13211`](https://github.com/noahlorinczcomi/simmrd/commit/8e13211e1055f7212343b1a7243c1b47f4835bac) [`fe52561`](https://github.com/noahlorinczcomi/simmrd/commit/fe52561b5b37abdc382d1ef667bf6ccf17cf3179) [`f75c6e0`](https://github.com/noahlorinczcomi/simmrd/commit/f75c6e0cd5343a78395f54f3d594296691b20a13))
- [`88fd437`](https://github.com/noahlorinczcomi/simmrd/commit/88fd437095b6f9b23795b228af83b6f38eedf923) Amend `nr`/`nc` to `nrow`/`ncol` to avoid partial argument matching
- [`30d971d`](https://github.com/noahlorinczcomi/simmrd/commit/30d971d3c5d5f7e1f816be54f00cdb598c39ba3d) Create .gitattributes
- [`c0f768c`](https://github.com/noahlorinczcomi/simmrd/commit/c0f768cc0dbc0442424d953d52f0eaf2ed185936) Convert CRLF line endings to LF
- [`4aa7c70`](https://github.com/noahlorinczcomi/simmrd/commit/4aa7c70d1170522a5235c3c9b3461a1c7fb034ed) Convert to LF line endings
- [`d50e5b5`](https://github.com/noahlorinczcomi/simmrd/commit/d50e5b503d2d65bb2d50043a54c617dd92b30e31) Add ggpubr to Imports
- [`b727d6b`](https://github.com/noahlorinczcomi/simmrd/commit/b727d6bdcc5a0ee511e3a0e28496f64700b06e21) Add ggplot2 and mvnfast to Imports
- [`251ebdd`](https://github.com/noahlorinczcomi/simmrd/commit/251ebdd29f0144121142d1beffd51212700b4c46) Build-ignore files in non-standard locations for an R package
- [`3cef46c`](https://github.com/noahlorinczcomi/simmrd/commit/3cef46cec48099ade6cff1af37a444075deb6d6b) usethis::use_mit_license()
- [`5cb96b0`](https://github.com/noahlorinczcomi/simmrd/commit/5cb96b037a0986b95ea3529158e7a8e5a5d91a39) Create .Rbuildignore
- [`ac2f540`](https://github.com/noahlorinczcomi/simmrd/commit/ac2f54027905e626af19726c4a8514318e25b244) Add Roxygen documentation generation project option
- [`700c50a`](https://github.com/noahlorinczcomi/simmrd/commit/700c50af078355b1904a5f5aca08d71612f4b043) usethis::use_tidy_description()
- [`b1cdd6d`](https://github.com/noahlorinczcomi/simmrd/commit/b1cdd6d02e0860775fc1aa2b0abdf9268b9367d5) Add BugReports entry
- [`8c6fa92`](https://github.com/noahlorinczcomi/simmrd/commit/8c6fa92bf0b706a1114c2abb0ae4f75e1a0c6f35) Add URL entry

## 2023-12-19

- [`dbbc78e`](https://github.com/noahlorinczcomi/simmrd/commit/dbbc78e3500aa048c28b20845c0bdc04f7cb7940) Update tutorial links in README from flat file to wiki
- [`7c88dda`](https://github.com/noahlorinczcomi/simmrd/commit/7c88ddab82831df56d0d077371652746b99da616) Add note to README that the paper changed substantially and is pending publication in Genetic Epidemiology

## 2023-12-11

- [`7212ea4`](https://github.com/noahlorinczcomi/simmrd/commit/7212ea477f00986e1735e6479e024208cc61160f) Comment out `true_variance_explained` output item (item 14) from README output list
- [`ac8fb36`](https://github.com/noahlorinczcomi/simmrd/commit/ac8fb361779c0cbfb3e74b26cdc71bda0b5e0fec) Add `simmr_flowchart_twodatagen.svg` (two-sample data generation flowchart)
- [`db97ffd`](https://github.com/noahlorinczcomi/simmrd/commit/db97ffd0d5c58df8e08846d87ce0084e9e3031fc) Replace the detailed `params` usage section in README with a Paper section placeholder
- [`7170f81`](https://github.com/noahlorinczcomi/simmrd/commit/7170f818964d6c375c198c34dea4653f89336b98) Remove obsolete content from `CHP_30kGWAS_100SNPs_1exposure_fullOverlap.R`
- [`dfe53d7`](https://github.com/noahlorinczcomi/simmrd/commit/dfe53d76820a9a49effe007d79dc9ef79146e070) Update `R/simmr.R`; add CHP 30kGWAS scenario script
- [`ea7dfa8`](https://github.com/noahlorinczcomi/simmrd/commit/ea7dfa835c81f6263f72352c4fa3d9838af29fef) Update README and tutorial; refresh example simulation and simtime figures

## 2023-12-06

- [`dfb8409`](https://github.com/noahlorinczcomi/simmrd/commit/dfb84094fd620eeb81d325a32cd85edbd25ce5b1) Add `example_params.R` (40-line example parameter configuration script)
- [`b5f08c1`](https://github.com/noahlorinczcomi/simmrd/commit/b5f08c1daf7f859ca74f125eee33f509944cd73f) Add `ld_figure.pdf`; clean up `params.R`
- [`3a83760`](https://github.com/noahlorinczcomi/simmrd/commit/3a83760e812b6c783eba7caf8f10995ea2a192ad) Add `example_sim_figure.pdf`/`.pptx`; update `plot_simdata_example` SVG

## 2023-12-05

- [`abb19b4`](https://github.com/noahlorinczcomi/simmrd/commit/abb19b46cb5e492bac98256d05b5b47213074918) Expand tutorial intro to explain `generate_individual()` vs `generate_summary()`
- [`5510ea6`](https://github.com/noahlorinczcomi/simmrd/commit/5510ea6bf6c7e3ff37ca06bcb9cd51580b5aaedf) Minor additions to tutorial.md
- [`0bb05d9`](https://github.com/noahlorinczcomi/simmrd/commit/0bb05d95fa8a62afd0f8b71d661aa9864798607f) Delete simmr directory
- [`fabe71a`](https://github.com/noahlorinczcomi/simmrd/commit/fabe71a51370de8747ad0ff7826e38d1ed56a485) Add `p1.svg` and `p2.svg` plot output files
- [`7fd30c4`](https://github.com/noahlorinczcomi/simmrd/commit/7fd30c4ea0feb7601c15e5ffd6e34e7629c1f5c8) Add initial `tutorial.md` (89 lines covering both `generate_summary()` and `generate_individual()`)
- [`f8491b4`](https://github.com/noahlorinczcomi/simmrd/commit/f8491b43ce27e967782edf3a821bcc6231448c9a) Add `R/simmr.R` (695 lines), man/ documentation stubs, and 80+ batch scenario scripts spanning CHP/UHP/WEAK × GWAS size × SNP count
- [`c6a9b8c`](https://github.com/noahlorinczcomi/simmrd/commit/c6a9b8caaeffee17ea4ef54a43e5a09da7f7991c) Add DESCRIPTION, NAMESPACE, `params.R`, timing CSVs, `plot_simdata_example` SVG, and `simmr.Rproj`
- [`5bdab60`](https://github.com/noahlorinczcomi/simmrd/commit/5bdab60ef677131b1554201ed42911d609499db7) Add full output list to README (items 1–14: `bx`, `by`, `bxse`, `byse`, `RhoME`, `LDMatrix`, etc.)
- [`a9a7248`](https://github.com/noahlorinczcomi/simmrd/commit/a9a7248b92e45ad57d64bc8fa515ce33f9234434) Remove old generate/params usage section from README; add Paper section placeholder

## 2023-11-22

- [`ca0cb77`](https://github.com/noahlorinczcomi/simmrd/commit/ca0cb7775b0367716fc99f60278a13400eb82231) Add `simmr/` package skeleton: R/ source, man/ pages, DESCRIPTION, NAMESPACE, LICENSE, and supporting scripts (`basicfunctions.R`, `generate_data.R`, `set_params.R`, etc.)

## 2023-09-13

- [`e930794`](https://github.com/noahlorinczcomi/simmrd/commit/e930794991863230cde3b4e27a5e15a91ea75768) Refactor `generate_data.R`
- [`8ee0b00`](https://github.com/noahlorinczcomi/simmrd/commit/8ee0b004525b5e93d1ba90ba558db5ee3fd18783) Refactor `basicfunctions.R`, `generate_data.R`, and `set_params.R`

## 2023-09-09

- [`0a2d9b8`](https://github.com/noahlorinczcomi/simmrd/commit/0a2d9b8c4193a96904ad4e2f50a9106c1604d814) Update `basicfunctions.R` and `set_params.R`
- [`5b8aedc`](https://github.com/noahlorinczcomi/simmrd/commit/5b8aedc62ca4a5ed87ab30e37083a3b8cb3fa872) Update `example_plotsimdata.svg` plot output
- [`1be502a`](https://github.com/noahlorinczcomi/simmrd/commit/1be502aa93897c7274691a5975eeb2af7f7ded0e) Add initial `example_plotsimdata.svg`
- [`ae9bd29`](https://github.com/noahlorinczcomi/simmrd/commit/ae9bd292aa2b75a3503a44d3f3aa34d0d201bfe1) Add helper functions to `basicfunctions.R`; refactor `generate_data.R` and `set_params.R`
- [`e84750d`](https://github.com/noahlorinczcomi/simmrd/commit/e84750d189b9d16675ea4656b7ab9e075fc1d3f0) Add `plot_simdata()` description and example SVG image to README

## 2023-09-08

- [`b3bcf0e`](https://github.com/noahlorinczcomi/simmrd/commit/b3bcf0ebe580c35fe84438a7503d695cb9753fa7) Update `generate_data.R`, `set_params.R`, and flowchart files

## 2023-09-01

- [`de951ac`](https://github.com/noahlorinczcomi/simmrd/commit/de951ac9d1beab3e475593a716196e081a72b27d) Update `generate_data.R`, `set_params.R`, and flowchart files

## 2023-08-31

- [`72cc4fc`](https://github.com/noahlorinczcomi/simmrd/commit/72cc4fc47ce4c7c4fa59b9545fd1e8bbc777622e) Update `simmr_flowchart` SVG and PDF
- [`255fa09`](https://github.com/noahlorinczcomi/simmrd/commit/255fa092589a6962236bd8155700090799374753) Add `.gitignore`; expand README tutorial with multi-method usage structure; add `signs_of_causal_effects` param; update `generate_data.R` and `set_params.R`
- [`d6b8794`](https://github.com/noahlorinczcomi/simmrd/commit/d6b8794c05aa06f44b24f96f1f64660a1a39aaec) Update `simmr_flowchart.svg`
- [`52d3199`](https://github.com/noahlorinczcomi/simmrd/commit/52d319978d9dc2eebf1c714a1ef16142c7b59d44) Add second usage example to README showing how to run simulations via `set_params.R`

## 2023-08-30

- [`225b06b`](https://github.com/noahlorinczcomi/simmrd/commit/225b06bdcfe5de956a2656b348dee492eb221230) Add detailed "Short tutorial" section to README with code example showing full parameter setup and usage

## 2023-08-29

- [`909f009`](https://github.com/noahlorinczcomi/simmrd/commit/909f0093748221d9af4b230cb1600e57d229e614) Update `simmr_flowchart` SVG and PDF
- [`5635fca`](https://github.com/noahlorinczcomi/simmrd/commit/5635fcfafef37d2dc527481da0f284ff1265b409) Add `simmr_flowchart.pdf`; update `simmr_flowchart.svg`

## 2023-08-28

- [`0450c87`](https://github.com/noahlorinczcomi/simmrd/commit/0450c879b6d2c74995ed9e3eee5176c26b4f3925) Update `simmr_flowchart.svg`

## 2023-08-27

- [`aa7ac9a`](https://github.com/noahlorinczcomi/simmrd/commit/aa7ac9ac2fcc79d32bfd054a9712410e533ddd8a) Clarify "scripts" → "R scripts" in README summary
- [`cc9f35b`](https://github.com/noahlorinczcomi/simmrd/commit/cc9f35b77fc2e8c9b8dfcec0fc3d05d985d56f88) Add flowchart SVG image to README
- [`d97cb37`](https://github.com/noahlorinczcomi/simmrd/commit/d97cb37079f13cd8ea21299153aabff6f13cc6f6) Rename `simmr_flowchartv4.svg` → `simmr_flowchart.svg`
- [`c7deec6`](https://github.com/noahlorinczcomi/simmrd/commit/c7deec61927a34a0508969efe25b2530ad336ad0) Add `simmr_flowchartv4.svg` (initial flowchart)
- [`59e9ea3`](https://github.com/noahlorinczcomi/simmrd/commit/59e9ea3deede04f2b117e23d7759c9c71963b15b) Remove `simmr-5.pdf`

## 2023-08-22

- [`64f0bfa`](https://github.com/noahlorinczcomi/simmrd/commit/64f0bfa2673407c318b5558d1cf1779bbf66a934) Fix line number reference for `source('generate_data.R')` in README (35 → 33)
- [`8d1e875`](https://github.com/noahlorinczcomi/simmrd/commit/8d1e8755609d5abdd16f9efbbfd8be12e36eb00e) Expand motivation text in README with citations to conflicting MR simulation papers; fix download link formatting
- [`1a3e80a`](https://github.com/noahlorinczcomi/simmrd/commit/1a3e80a0b8c9e3330f631945eb1c6ec3166a5a0c) Simplify `basicfunctions.R`, `generate_data.R`, and `set_params.R` (cleanup)

## 2023-08-21

- [`2afd90d`](https://github.com/noahlorinczcomi/simmrd/commit/2afd90dda72b1dd4980efb2f9608490eeebab375) Add `simmr-5.pdf` (working paper draft — not peer-checked or reviewed)
- [`25188d3`](https://github.com/noahlorinczcomi/simmrd/commit/25188d3fdd1a004488069a0dc40f68ee04745253) Fix `generate_data.R`: set genetic effect mean to 0 and add genetic correlation check via Cholesky decomposition

## 2023-08-20

- [`b70f26a`](https://github.com/noahlorinczcomi/simmrd/commit/b70f26afed0d7a87ccc558aabd4b2ffc917a7301) Update all core R scripts: `basicfunctions.R`, `creating_plots.R`, `generate_data.R`, `set_params.R`, `timing_datagen.R`
- [`13d0c95`](https://github.com/noahlorinczcomi/simmrd/commit/13d0c954d2b7bfb18ed416f79834edb6d8bec98e) Update README output list: add `IVtype`, `bxunstd`, `bxseunstd`; replace `mSelected`/`mNotPruned` with `mIVs`

## 2023-08-17

- [`214eec0`](https://github.com/noahlorinczcomi/simmrd/commit/214eec04257fbc90d6c5be465da1a1ece7a42fc4) Add functions to `basicfunctions.R` and `creating_plots.R`; update `generate_data.R` and `set_params.R`

## 2023-08-16

- [`e6a65ff`](https://github.com/noahlorinczcomi/simmrd/commit/e6a65ffc7567fab83b0366f3f50eb3e0b5e0add1) Remove `fig.pdf` and `fig.pptx`
- [`3d34179`](https://github.com/noahlorinczcomi/simmrd/commit/3d3417921e2e5e7791355026b94d85e3d58e9d5c) Add `fig.pdf`; update `fig.pptx`
- [`667281b`](https://github.com/noahlorinczcomi/simmrd/commit/667281bb9181be471b004eddfea636031b2cc28b) Add `fig.pptx`; update `generate_data.R`
- [`85a71b1`](https://github.com/noahlorinczcomi/simmrd/commit/85a71b1fbf7d23aed1d517a2b991b9ea20c32a0c) Update `generate_data.R` and `set_params.R`
- [`a51a080`](https://github.com/noahlorinczcomi/simmrd/commit/a51a0807a99e286d4463a4ac1a81ef3ffecab41f) Update `generate_data.R` and `set_params.R`
- [`707ba10`](https://github.com/noahlorinczcomi/simmrd/commit/707ba103c061f9efc5de970260d3ce200be935e0) Update `set_params.R` parameters; add `timing_datagen.R`

## 2023-08-15

- [`9feb91b`](https://github.com/noahlorinczcomi/simmrd/commit/9feb91b9d673fa9e3dd98a68bbc1c931d0cd7253) Create README.md
- [`4c91a4d`](https://github.com/noahlorinczcomi/simmrd/commit/4c91a4dbcf4f40d3d71e03f3acb2bb403ddbb7aa) Add initial core scripts: `basicfunctions.R`, `generate_data.R`, `set_params.R`, `simdata.Rds`, and `simmr.Rproj`
- [`b92114c`](https://github.com/noahlorinczcomi/simmrd/commit/b92114c4d1ef3699f7f829e26af2bdb16a441926) Delete `simdata.Rds`
- [`17987362`](https://github.com/noahlorinczcomi/simmrd/commit/17987362a5f2638e197760d39e5e4284afa8f89c) Large revisions to `basicfunctions.R`, `generate_data.R`, and `set_params.R`
- [`b4bbe4d`](https://github.com/noahlorinczcomi/simmrd/commit/b4bbe4d06105c37db4333af0b792bf267e98b636) Update `generate_data.R` and `set_params.R`; add `simmr.Rproj`
- [`2e51642`](https://github.com/noahlorinczcomi/simmrd/commit/2e5164244d7470d603b00ed71de684647dad792b) Update `set_params.R`
- [`fbb88ec`](https://github.com/noahlorinczcomi/simmrd/commit/fbb88eca291d5be95edc578038eeb0481c267c82) Expand "Who is this for?" section; add simulation scenario list (sample overlap, UHP, CHP, weak instruments, winner's curse, LD) to README
- [`f435721`](https://github.com/noahlorinczcomi/simmrd/commit/f4357215986ae194ecf2afe2e5a83ec1cc794c0d) Fix reference link formatting; add full "How do I use it?" walkthrough with output object descriptions
- [`faeb741`](https://github.com/noahlorinczcomi/simmrd/commit/faeb741e275e68163858e8e910f01a0b1b8175d0) Clarify output object descriptions in README; reformat "Can I add?" section
- [`4d5a4aa`](https://github.com/noahlorinczcomi/simmrd/commit/4d5a4aa52085692aface639d5693b3a3b1482753) Fix typo in README open-source notice ("but" → "by")
- [`a29413e`](https://github.com/noahlorinczcomi/simmrd/commit/a29413ec6f4d670c6fc9bc555028a1c0f6d6a052) Clarify working directory instruction in README usage steps
- [`1a7bbee`](https://github.com/noahlorinczcomi/simmrd/commit/1a7bbeed32b92fe3208dcf64c3bdf220abe4b348) Add `theta` to output list in README; fix output object count (9 → 10)
- [`5941cf3`](https://github.com/noahlorinczcomi/simmrd/commit/5941cf30b4f4b8ac2cb1e2f09df0c319a113879b) Initial commit
