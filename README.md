<p align="center">
<img src=https://github.com/MichelNivard/trigenometry/assets/11858442/63d46a1c-39d6-47d5-b49c-082f0c960874>
</p>

Non linear genetic dependence using LD score regression and high school trigonometry


### Installation

Until we go public contributors can make a github auth token:

https://github.com/settings/tokens

And then install from github:

```r
devtools::install_github("WonuAkingbuwa/trigenometry",
                         ref = "main",
                         auth_token = # YOUR TOKEN HERE, DO NOT EVER SHARE TOKEN                      
                        )

```


### Usage

Currently we split out making a curve from plotting it:

```r
c2c <- cor2curve(...)

plot(c2c)

```

This currently allows the user to set some axis labels, It will ease future updates for manually setting colors and perhaps ease combining plots?


