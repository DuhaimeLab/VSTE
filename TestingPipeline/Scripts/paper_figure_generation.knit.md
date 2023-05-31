
<!-- rnb-text-begin -->

---
title: "Pipeline Performance Figure Generation"
output: html_notebook
---

This R script creates all of the figures used to compare the performance across pipelines.

This Rmarkdown file assesses the output of CheckV, DeepVirFinder, Kaiju,
VIBRANT, VirSorter, and VirSorter2 on multiple training sets of microbial DNA, 
primarily from NCBI. Created from fungal, viral, bacterial, archeael, protist,
and plasmid DNA sequences

Please reach out to James Riddell (riddell.26@buckeyemail.osu.edu) or
Bridget Hegarty (beh53@case.edu) regarding any issues, or open an issue on github.


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShnZ3Bsb3QyKVxuYGBgIn0= -->

```r
library(ggplot2)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiVGhlcmUgd2VyZSA1MCBvciBtb3JlIHdhcm5pbmdzICh1c2Ugd2FybmluZ3MoKSB0byBzZWUgdGhlIGZpcnN0IDUwKVxuIn0= -->

```
There were 50 or more warnings (use warnings() to see the first 50)
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShwbHlyKVxubGlicmFyeShyZXNoYXBlMilcbmxpYnJhcnkodmlyaWRpcylcbmxpYnJhcnkodGlkeXIpXG5saWJyYXJ5KGRwbHlyKVxubGlicmFyeShyZWFkcilcbmxpYnJhcnkoZGF0YS50YWJsZSlcbmxpYnJhcnkocFJPQylcbmxpYnJhcnkoXCJzdHJpbmdyXCIpXG5gYGAifQ== -->

```r
library(plyr)
library(reshape2)
library(viridis)
library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(pROC)
library("stringr")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Import the file that combines the results from each of the tools from running "combining_tool_output.Rmd":

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyA8LSByZWFkX3RzdihcIi4uL0ludGVybWVkaWFyeUZpbGVzL3ZpcmFsX3Rvb2xzX2NvbWJpbmVkLnRzdlwiKVxuYGBgIn0= -->

```r
viruses <- read_tsv("../IntermediaryFiles/viral_tools_combined.tsv")
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiXG7ilIDilIAgQ29sdW1uIHNwZWNpZmljYXRpb24g4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSA4pSAXG5jb2xzKFxuICAuZGVmYXVsdCA9IGNvbF9kb3VibGUoKSxcbiAgc2VxdHlwZSA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgY29udGlnID0gY29sX2NoYXJhY3RlcigpLFxuICBjaGVja3ZfcHJvdmlydXMgPSBjb2xfY2hhcmFjdGVyKCksXG4gIGNoZWNrdl9xdWFsaXR5ID0gY29sX2NoYXJhY3RlcigpLFxuICBtZXRob2QueCA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgQ2xhc3NpZmllZCA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgSURzX2FsbCA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgU2VxID0gY29sX2NoYXJhY3RlcigpLFxuICBLYWlqdV9WaXJhbCA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgS2luZ2RvbSA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgdHlwZSA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgdmlicmFudF9xdWFsaXR5ID0gY29sX2NoYXJhY3RlcigpLFxuICBtZXRob2QueSA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgdmlicmFudF9wcm9waGFnZSA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgdnMydHlwZSA9IGNvbF9jaGFyYWN0ZXIoKSxcbiAgbWF4X3Njb3JlX2dyb3VwID0gY29sX2NoYXJhY3RlcigpLFxuICBwcm92aXJ1cyA9IGNvbF9sb2dpY2FsKClcbilcbuKEuSBVc2UgYHNwZWMoKWAgZm9yIHRoZSBmdWxsIGNvbHVtbiBzcGVjaWZpY2F0aW9ucy5cbiJ9 -->

```

── Column specification ────────────────────────────────────────────────────────────────────────────────────
cols(
  .default = col_double(),
  seqtype = col_character(),
  contig = col_character(),
  checkv_provirus = col_character(),
  checkv_quality = col_character(),
  method.x = col_character(),
  Classified = col_character(),
  IDs_all = col_character(),
  Seq = col_character(),
  Kaiju_Viral = col_character(),
  Kingdom = col_character(),
  type = col_character(),
  vibrant_quality = col_character(),
  method.y = col_character(),
  vibrant_prophage = col_character(),
  vs2type = col_character(),
  max_score_group = col_character(),
  provirus = col_logical()
)
ℹ Use `spec()` for the full column specifications.
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrYWlqdV9tYXRjaF9yYXRpbyA8LSB2aXJ1c2VzJGxlbi92aXJ1c2VzJGNoZWNrdl9sZW5ndGhcbmBgYCJ9 -->

```r
viruses$kaiju_match_ratio <- viruses$len/viruses$checkv_length
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



Function that allows for the comparing pieces of the pipeline

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2V0dGluZ192aXJhbF9zZXRfMSA8LSBmdW5jdGlvbihpbnB1dF9zZXFzLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQ9RkFMU0UsIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjI9RkFMU0UsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfZGVlcHZpcmZpbmRlcj1GQUxTRSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWw9RkFMU0UsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR2XzE9VCwgdHZfMj1ULCB0dl8zPVQsIHR2XzQ9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsPUZBTFNFLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBudHZfMT1ULCBudHZfMj1ULCBudHZfMz1ULCBudHZfND1ULCBudHZfNT1GLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcj1GQUxTRSkge1xuICBcbiAga2VlcF9zY29yZSA8LSByZXAoMCwgbnJvdyhpbnB1dF9zZXFzKSlcblxuICBpZiAoaW5jbHVkZV92aWJyYW50KSB7XG4gICAga2VlcF9zY29yZVtpbnB1dF9zZXFzJHZpYnJhbnRfcXVhbGl0eT09XCJoaWdoIHF1YWxpdHkgZHJhZnRcIl0gPC0ga2VlcF9zY29yZVtpbnB1dF9zZXFzJHZpYnJhbnRfcXVhbGl0eT09XCJoaWdoIHF1YWxpdHkgZHJhZnRcIl0gKyAxXG4gICAga2VlcF9zY29yZVtpbnB1dF9zZXFzJHZpYnJhbnRfcXVhbGl0eT09XCJtZWRpdW0gcXVhbGl0eSBkcmFmdFwiXSA8LSBrZWVwX3Njb3JlW2lucHV0X3NlcXMkdmlicmFudF9xdWFsaXR5PT1cIm1lZGl1bSBxdWFsaXR5IGRyYWZ0XCJdICsgMVxuICAgIGtlZXBfc2NvcmVbaW5wdXRfc2VxcyR2aWJyYW50X3F1YWxpdHk9PVwibG93IHF1YWxpdHkgZHJhZnRcIiAmIGlucHV0X3NlcXMkcHJvdmlydXNdIDwtIGtlZXBfc2NvcmVbaW5wdXRfc2VxcyR2aWJyYW50X3F1YWxpdHk9PVwibG93IHF1YWxpdHkgZHJhZnRcIiAmIGlucHV0X3NlcXMkcHJvdmlydXNdICsgMC41XG4gIH1cbiAgXG4gIGlmIChpbmNsdWRlX3ZpcnNvcnRlcjIpIHtcbiAgICBrZWVwX3Njb3JlW2lucHV0X3NlcXMkbWF4X3Njb3JlPj0wLjUwXSA8LSBrZWVwX3Njb3JlW2lucHV0X3NlcXMkbWF4X3Njb3JlPj0wLjUwXSArIDAuNVxuICAgIGtlZXBfc2NvcmVbaW5wdXRfc2VxcyRtYXhfc2NvcmU+PTAuOTVdIDwtIGtlZXBfc2NvcmVbaW5wdXRfc2VxcyRtYXhfc2NvcmU+PTAuOTVdICsgMC41ICBcbiAgfVxuICBcbiAgaWYgKGluY2x1ZGVfdmlyc29ydGVyKSB7XG4gICAga2VlcF9zY29yZVtpbnB1dF9zZXFzJGNhdGVnb3J5PT0xXSA8LSBrZWVwX3Njb3JlW2lucHV0X3NlcXMkY2F0ZWdvcnk9PTFdICsgMVxuICAgIGtlZXBfc2NvcmVbaW5wdXRfc2VxcyRjYXRlZ29yeT09Ml0gPC0ga2VlcF9zY29yZVtpbnB1dF9zZXFzJGNhdGVnb3J5PT0yXSArIDAuNVxuICAgIGtlZXBfc2NvcmVbaW5wdXRfc2VxcyRjYXRlZ29yeT09M10gPC0ga2VlcF9zY29yZVtpbnB1dF9zZXFzJGNhdGVnb3J5PT0zXSArIDAuNVxuICAgIGtlZXBfc2NvcmVbaW5wdXRfc2VxcyRjYXRlZ29yeT09NF0gPC0ga2VlcF9zY29yZVtpbnB1dF9zZXFzJGNhdGVnb3J5PT00XSArIDFcbiAgICBrZWVwX3Njb3JlW2lucHV0X3NlcXMkY2F0ZWdvcnk9PTVdIDwtIGtlZXBfc2NvcmVbaW5wdXRfc2VxcyRjYXRlZ29yeT09NV0gKyAwLjUgXG4gICAga2VlcF9zY29yZVtpbnB1dF9zZXFzJGNhdGVnb3J5PT02XSA8LSBrZWVwX3Njb3JlW2lucHV0X3NlcXMkY2F0ZWdvcnk9PTZdICsgMC41IFxuICB9XG4gIFxuICBpZiAoaW5jbHVkZV9kZWVwdmlyZmluZGVyKSB7XG4gICAgICMgYWRkIGlmIERWRiBjYWxscyB2aXJhbFxuICAgIGtlZXBfc2NvcmVbKGlucHV0X3NlcXMkc2NvcmU+PTAuOSAmIGlucHV0X3NlcXMkcHZhbHVlPD0wLjA1KSAmIGlucHV0X3NlcXMkY2hlY2t2X2xlbmd0aDwyMDAwMF0gPC0ga2VlcF9zY29yZVsoaW5wdXRfc2VxcyRzY29yZT49MC45ICYgaW5wdXRfc2VxcyRwdmFsdWU8PTAuMDUpICYgaW5wdXRfc2VxcyRjaGVja3ZfbGVuZ3RoPDIwMDAwXSArIDAuNVxuICAgIGtlZXBfc2NvcmVbKGlucHV0X3NlcXMkc2NvcmU+PTAuNyAmIGlucHV0X3NlcXMkcHZhbHVlPD0wLjA1KSAmIGlucHV0X3NlcXMkY2hlY2t2X2xlbmd0aDwyMDAwMF0gPC0ga2VlcF9zY29yZVsoaW5wdXRfc2VxcyRzY29yZT49MC43ICYgaW5wdXRfc2VxcyRwdmFsdWU8PTAuMDUpICYgaW5wdXRfc2VxcyRjaGVja3ZfbGVuZ3RoPDIwMDAwXSArIDAuNVxuICB9XG4gIFxuICBpZiAoaW5jbHVkZV90dW5pbmdfdmlyYWwpIHtcbiAgICBrZWVwX3Njb3JlIDwtIHR1bmluZ192aXJhbChjdXJyZW50X3Njb3JlPWtlZXBfc2NvcmUsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5zPWlucHV0X3NlcXMsXG4gICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dl8xPXR2XzEsXG4gICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dl8yPXR2XzIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dl8zPXR2XzMsXG4gICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dl80PXR2XzQpXG4gIH1cbiAgXG4gIGlmIChpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwpIHtcbiAgICAjIHR1bmluZyByZW1vdmFsIFxuICAgIGtlZXBfc2NvcmUgPC0gdHVuaW5nX25vdF92aXJhbChjdXJyZW50X3Njb3JlPWtlZXBfc2NvcmUsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5zPWlucHV0X3NlcXMsXG4gICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dl8xPW50dl8xLFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfMj1udHZfMixcbiAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R2XzM9bnR2XzMsXG4gICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dl80PW50dl80LFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfNT1udHZfNSlcbiAgIFxuICB9XG4gIFxuICByZXR1cm4oa2VlcF9zY29yZSlcbiAgXG59XG5cbnR1bmluZ192aXJhbCA8LSBmdW5jdGlvbihjdXJyZW50X3Njb3JlLCBpbnMsIFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfMT1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfMj1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfMz1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfND1UKSB7XG4gICAgIyB0dW5pbmcgYWRkaXRpb25cbiAgaWYgKGluY2x1ZGVfdHZfMSkge1xuICAgIGN1cnJlbnRfc2NvcmVbaW5zJGhhbGxtYXJrPjJdIDwtIGN1cnJlbnRfc2NvcmVbaW5zJGhhbGxtYXJrPjJdICsgMVxuICB9XG4gIGlmIChpbmNsdWRlX3R2XzIpIHtcbiAgICBjdXJyZW50X3Njb3JlW2lucyRLYWlqdV9WaXJhbD09XCJWaXJ1c2VzXCIgJiBpbnMka2FpanVfbWF0Y2hfcmF0aW8+PTAuM10gPC0gY3VycmVudF9zY29yZVtpbnMkS2FpanVfVmlyYWw9PVwiVmlydXNlc1wiICYgaW5zJGthaWp1X21hdGNoX3JhdGlvPj0wLjNdICsgMVxuICB9XG4gICAgI25vdGU6IGhhZCB0cmllZCBwdWxsaW5nIG91dCB0aGlzIHJ1bGUsIGJ1dCBiZXR0ZXIgcmVjYWxsIHdpdGhvdXQgc2FjcmlmaWNpbmcgcHJlY2lzaW9uIHdpdGggaXQgaW5cbiAgICAjbm90ZTogaGF2aW5nIGthaWp1IHJlbW92YWwgcnVsZSBkaWRuJ3QgaGVscFxuICBpZiAoaW5jbHVkZV90dl8zKSB7XG4gICAgY3VycmVudF9zY29yZVtpbnMkcGVyY2VudF91bmtub3duPj03NSAmIGlucyRjaGVja3ZfbGVuZ3RoPDUwMDAwXSA8LSBjdXJyZW50X3Njb3JlW2lucyRwZXJjZW50X3Vua25vd24+PTc1ICYgaW5zJGNoZWNrdl9sZW5ndGg8NTAwMDBdICsgMC41XG4gIH1cbiAgaWYgKGluY2x1ZGVfdHZfNCkge1xuICAgIGN1cnJlbnRfc2NvcmVbaW5zJHZpcmFsPj01MCB8IGlucyRwZXJjZW50X3ZpcmFsPj01MF0gPC0gY3VycmVudF9zY29yZVtpbnMkdmlyYWw+PTUwIHwgaW5zJHBlcmNlbnRfdmlyYWw+PTUwXSArIDAuNSBcbiAgfVxuICByZXR1cm4oY3VycmVudF9zY29yZSlcbn1cblxudHVuaW5nX25vdF92aXJhbCA8LSBmdW5jdGlvbihjdXJyZW50X3Njb3JlLCBpbnMsIFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfMT1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfMj1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfMz1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfND1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHZfNT1GKSB7XG4gICAgIyB0dW5pbmcgcmVtb3ZhbFxuICBpZiAoaW5jbHVkZV90dl8yKSB7XG4gICAgY3VycmVudF9zY29yZVtpbnMkY2hlY2t2X3ZpcmFsX2dlbmVzPT0wICYgaW5zJGNoZWNrdl9ob3N0X2dlbmVzPj0xXSA8LSBjdXJyZW50X3Njb3JlW2lucyRjaGVja3ZfdmlyYWxfZ2VuZXM9PTAgJiBpbnMkY2hlY2t2X2hvc3RfZ2VuZXM+PTFdIC0gMVxuICB9XG4gIGlmIChpbmNsdWRlX3R2XzMpIHtcbiAgICBjdXJyZW50X3Njb3JlWygoaW5zJGNoZWNrdl92aXJhbF9nZW5lcyozKSA8PSBpbnMkY2hlY2t2X2hvc3RfZ2VuZXMpICYgIWlucyRwcm92aXJ1c10gPC0gY3VycmVudF9zY29yZVsoKGlucyRjaGVja3ZfdmlyYWxfZ2VuZXMqMykgPD0gaW5zJGNoZWNrdl9ob3N0X2dlbmVzKSAmICFpbnMkcHJvdmlydXNdIC0gMVxuICB9XG4gIFxuICBpZiAoaW5jbHVkZV90dl81KSB7XG4gICAgIyByZW1vdmUgaWYgRFZGIGNhbGxzIGl0IG5vdCB2aXJhbFxuICAgIGN1cnJlbnRfc2NvcmVbaW5zJHNjb3JlPD0wLjcgJiBpbnMkcHZhbHVlPD0wLjA1XSA8LSBjdXJyZW50X3Njb3JlW2lucyRzY29yZTw9MC43ICYgaW5zJHB2YWx1ZTw9MC4wNV0gLSAxIFxuICB9XG4gIFxuICBpZiAoaW5jbHVkZV90dl8xKSB7XG4gICAgY3VycmVudF9zY29yZVsoaW5zJGNoZWNrdl9ob3N0X2dlbmVzPjUwKSAmICFpbnMkcHJvdmlydXNdIDwtIC0zXG4gIH1cbiAgXG4gIGlmIChpbmNsdWRlX3R2XzQpIHtcbiAgICBjdXJyZW50X3Njb3JlW2lucyRjaGVja3ZfbGVuZ3RoPjUwMDAwMCAmIGlucyRoYWxsbWFyazw9MV0gPC0gLTNcbiAgfVxuICBcbiAgcmV0dXJuKGN1cnJlbnRfc2NvcmUpXG59XG5gYGAifQ== -->

```r
getting_viral_set_1 <- function(input_seqs,
                                include_vibrant=FALSE, 
                                include_virsorter2=FALSE,
                                include_deepvirfinder=FALSE,
                                include_tuning_viral=FALSE,
                                tv_1=T, tv_2=T, tv_3=T, tv_4=T,
                                include_tuning_not_viral=FALSE,
                                ntv_1=T, ntv_2=T, ntv_3=T, ntv_4=T, ntv_5=F,
                                include_virsorter=FALSE) {
  
  keep_score <- rep(0, nrow(input_seqs))

  if (include_vibrant) {
    keep_score[input_seqs$vibrant_quality=="high quality draft"] <- keep_score[input_seqs$vibrant_quality=="high quality draft"] + 1
    keep_score[input_seqs$vibrant_quality=="medium quality draft"] <- keep_score[input_seqs$vibrant_quality=="medium quality draft"] + 1
    keep_score[input_seqs$vibrant_quality=="low quality draft" & input_seqs$provirus] <- keep_score[input_seqs$vibrant_quality=="low quality draft" & input_seqs$provirus] + 0.5
  }
  
  if (include_virsorter2) {
    keep_score[input_seqs$max_score>=0.50] <- keep_score[input_seqs$max_score>=0.50] + 0.5
    keep_score[input_seqs$max_score>=0.95] <- keep_score[input_seqs$max_score>=0.95] + 0.5  
  }
  
  if (include_virsorter) {
    keep_score[input_seqs$category==1] <- keep_score[input_seqs$category==1] + 1
    keep_score[input_seqs$category==2] <- keep_score[input_seqs$category==2] + 0.5
    keep_score[input_seqs$category==3] <- keep_score[input_seqs$category==3] + 0.5
    keep_score[input_seqs$category==4] <- keep_score[input_seqs$category==4] + 1
    keep_score[input_seqs$category==5] <- keep_score[input_seqs$category==5] + 0.5 
    keep_score[input_seqs$category==6] <- keep_score[input_seqs$category==6] + 0.5 
  }
  
  if (include_deepvirfinder) {
     # add if DVF calls viral
    keep_score[(input_seqs$score>=0.9 & input_seqs$pvalue<=0.05) & input_seqs$checkv_length<20000] <- keep_score[(input_seqs$score>=0.9 & input_seqs$pvalue<=0.05) & input_seqs$checkv_length<20000] + 0.5
    keep_score[(input_seqs$score>=0.7 & input_seqs$pvalue<=0.05) & input_seqs$checkv_length<20000] <- keep_score[(input_seqs$score>=0.7 & input_seqs$pvalue<=0.05) & input_seqs$checkv_length<20000] + 0.5
  }
  
  if (include_tuning_viral) {
    keep_score <- tuning_viral(current_score=keep_score,
                               ins=input_seqs,
                         include_tv_1=tv_1,
                         include_tv_2=tv_2,
                         include_tv_3=tv_3,
                         include_tv_4=tv_4)
  }
  
  if (include_tuning_not_viral) {
    # tuning removal 
    keep_score <- tuning_not_viral(current_score=keep_score,
                               ins=input_seqs,
                         include_tv_1=ntv_1,
                         include_tv_2=ntv_2,
                         include_tv_3=ntv_3,
                         include_tv_4=ntv_4,
                         include_tv_5=ntv_5)
   
  }
  
  return(keep_score)
  
}

tuning_viral <- function(current_score, ins, 
                         include_tv_1=T,
                         include_tv_2=T,
                         include_tv_3=T,
                         include_tv_4=T) {
    # tuning addition
  if (include_tv_1) {
    current_score[ins$hallmark>2] <- current_score[ins$hallmark>2] + 1
  }
  if (include_tv_2) {
    current_score[ins$Kaiju_Viral=="Viruses" & ins$kaiju_match_ratio>=0.3] <- current_score[ins$Kaiju_Viral=="Viruses" & ins$kaiju_match_ratio>=0.3] + 1
  }
    #note: had tried pulling out this rule, but better recall without sacrificing precision with it in
    #note: having kaiju removal rule didn't help
  if (include_tv_3) {
    current_score[ins$percent_unknown>=75 & ins$checkv_length<50000] <- current_score[ins$percent_unknown>=75 & ins$checkv_length<50000] + 0.5
  }
  if (include_tv_4) {
    current_score[ins$viral>=50 | ins$percent_viral>=50] <- current_score[ins$viral>=50 | ins$percent_viral>=50] + 0.5 
  }
  return(current_score)
}

tuning_not_viral <- function(current_score, ins, 
                         include_tv_1=T,
                         include_tv_2=T,
                         include_tv_3=T,
                         include_tv_4=T,
                         include_tv_5=F) {
    # tuning removal
  if (include_tv_2) {
    current_score[ins$checkv_viral_genes==0 & ins$checkv_host_genes>=1] <- current_score[ins$checkv_viral_genes==0 & ins$checkv_host_genes>=1] - 1
  }
  if (include_tv_3) {
    current_score[((ins$checkv_viral_genes*3) <= ins$checkv_host_genes) & !ins$provirus] <- current_score[((ins$checkv_viral_genes*3) <= ins$checkv_host_genes) & !ins$provirus] - 1
  }
  
  if (include_tv_5) {
    # remove if DVF calls it not viral
    current_score[ins$score<=0.7 & ins$pvalue<=0.05] <- current_score[ins$score<=0.7 & ins$pvalue<=0.05] - 1 
  }
  
  if (include_tv_1) {
    current_score[(ins$checkv_host_genes>50) & !ins$provirus] <- -3
  }
  
  if (include_tv_4) {
    current_score[ins$checkv_length>500000 & ins$hallmark<=1] <- -3
  }
  
  return(current_score)
}
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


# Assessing performance against the "truth"
note that this is only as accurate as the annotations of the input sequences

this function calculates the precision, recall, and F1 score for each pipeline

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYXNzZXNzX3BlcmZvcm1hbmNlIDwtIGZ1bmN0aW9uKHNlcXR5cGUsIGtlZXBfc2NvcmUpIHtcbiAgXG4gIHRydWVwb3NpdGl2ZSA8LSByZXAoXCJub3QgdmlyYWxcIiwgbGVuZ3RoKHNlcXR5cGUpKVxuICB0cnVlcG9zaXRpdmVbc2VxdHlwZT09XCJ2aXJ1c1wiXSA8LSBcInZpcmFsXCJcbiAgXG4gICNtYWtlIGNvbmZ1c2lvbiBtYXRyaXhcbiAgY29uZnVzaW9uX21hdHJpeCA8LSByZXAoXCJ0cnVlIG5lZ2F0aXZlXCIsIGxlbmd0aChrZWVwX3Njb3JlKSlcbiAgY29uZnVzaW9uX21hdHJpeFt0cnVlcG9zaXRpdmU9PVwidmlyYWxcIiAmIGtlZXBfc2NvcmU8PTFdIDwtIFwiZmFsc2UgbmVnYXRpdmVcIlxuICBjb25mdXNpb25fbWF0cml4W3RydWVwb3NpdGl2ZT09XCJ2aXJhbFwiICYga2VlcF9zY29yZT49MV0gPC0gXCJ0cnVlIHBvc2l0aXZlXCJcbiAgY29uZnVzaW9uX21hdHJpeFt0cnVlcG9zaXRpdmU9PVwibm90IHZpcmFsXCIgJiBrZWVwX3Njb3JlPj0xXSA8LSBcImZhbHNlIHBvc2l0aXZlXCJcbiAgXG4gIFRQIDwtIHRhYmxlKGNvbmZ1c2lvbl9tYXRyaXgpWzRdXG4gIEZQIDwtIHRhYmxlKGNvbmZ1c2lvbl9tYXRyaXgpWzJdXG4gIFROIDwtIHRhYmxlKGNvbmZ1c2lvbl9tYXRyaXgpWzNdXG4gIEZOIDwtIHRhYmxlKGNvbmZ1c2lvbl9tYXRyaXgpWzFdXG4gIFxuICBwcmVjaXNpb24gPC0gVFAvKFRQK0ZQKVxuICByZWNhbGwgPC0gVFAvKFRQK0ZOKVxuICBGMSA8LSAyKnByZWNpc2lvbipyZWNhbGwvKHByZWNpc2lvbityZWNhbGwpXG4gIFxuICBNQ0MgPC0gKFRQKlROLUZQKkZOKS9zcXJ0KGFzLm51bWVyaWMoVFArRlApKmFzLm51bWVyaWMoVFArRk4pKmFzLm51bWVyaWMoVE4rRlApKmFzLm51bWVyaWMoVE4rRk4pKVxuICBcbiAgcHJvcF92aXJhbCA8LSAoVFArRlApLyhUUCtGUCtUTitGTilcbiAgXG4gIHBlcmZvcm1hbmNlIDwtIGMocHJlY2lzaW9uLCByZWNhbGwsIEYxLCBNQ0MsIHByb3BfdmlyYWwpXG4gIG5hbWVzKHBlcmZvcm1hbmNlKSA8LSBjKFwicHJlY2lzaW9uXCIsIFwicmVjYWxsXCIsIFwiRjFcIiwgXCJNQ0NcIiwgXCJwcm9wX3ZpcmFsXCIpXG4gIFxuICByZXR1cm4ocGVyZm9ybWFuY2UpXG59XG5gYGAifQ== -->

```r
assess_performance <- function(seqtype, keep_score) {
  
  truepositive <- rep("not viral", length(seqtype))
  truepositive[seqtype=="virus"] <- "viral"
  
  #make confusion matrix
  confusion_matrix <- rep("true negative", length(keep_score))
  confusion_matrix[truepositive=="viral" & keep_score<=1] <- "false negative"
  confusion_matrix[truepositive=="viral" & keep_score>=1] <- "true positive"
  confusion_matrix[truepositive=="not viral" & keep_score>=1] <- "false positive"
  
  TP <- table(confusion_matrix)[4]
  FP <- table(confusion_matrix)[2]
  TN <- table(confusion_matrix)[3]
  FN <- table(confusion_matrix)[1]
  
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- 2*precision*recall/(precision+recall)
  
  MCC <- (TP*TN-FP*FN)/sqrt(as.numeric(TP+FP)*as.numeric(TP+FN)*as.numeric(TN+FP)*as.numeric(TN+FN))
  
  prop_viral <- (TP+FP)/(TP+FP+TN+FN)
  
  performance <- c(precision, recall, F1, MCC, prop_viral)
  names(performance) <- c("precision", "recall", "F1", "MCC", "prop_viral")
  
  return(performance)
}
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


combination of tools list

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY29tYm9zX2xpc3QgPC0gZGF0YS5mcmFtZSh0b29sY29tYm89cmVwKDAsIDY0KSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgdHVuZV9ub3RfdmlyYWw9cmVwKDAsIDY0KSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgRFZGPXJlcCgwLCA2NCksXG4gICAgICAgICAgICAgICAgICAgICAgICAgIHR1bmVfdmlyYWw9cmVwKDAsIDY0KSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgVklCUkFOVD1yZXAoMCwgNjQpLFxuICAgICAgICAgICAgICAgICAgICAgICAgICBWUz1yZXAoMCwgNjQpLFxuICAgICAgICAgICAgICAgICAgICAgICAgICBWUzI9cmVwKDAsIDY0KSlcbnAgPC0gMVxuXG5mb3IgKGkgaW4gYygwLDEpKXtcbiAgZm9yIChqIGluIGMoMCwxKSl7XG4gICAgZm9yIChrIGluIGMoMCwxKSl7XG4gICAgICBmb3IgKGwgaW4gYygwLDEpKXtcbiAgICAgICAgZm9yIChtIGluIGMoMCwxKSl7XG4gICAgICAgICAgZm9yIChuIGluIGMoMCwxKSl7XG4gICAgICAgICAgICBjb21ib3NfbGlzdCR0b29sY29tYm9bcF0gPC0gcGFzdGUoaSxqLGssbCxtLG4pXG4gICAgICAgICAgICBjb21ib3NfbGlzdCR0b29sY29tYm8yW3BdIDwtIHBhc3RlKGlmKGkpe1widHZcIn1lbHNle1wiMFwifSxpZihqKXtcIkRWRlwifWVsc2V7XCIwXCJ9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZihrKXtcInRudlwifWVsc2V7XCIwXCJ9LGlmKGwpe1wiVkJcIn1lbHNle1wiMFwifSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaWYobSl7XCJWU1wifWVsc2V7XCIwXCJ9LGlmKG4pe1wiVlMyXCJ9ZWxzZXtcIjBcIn0pXG4gICAgICAgICAgICBjb21ib3NfbGlzdCR0dW5lX25vdF92aXJhbFtwXSA8LSBrXG4gICAgICAgICAgICBjb21ib3NfbGlzdCREVkZbcF0gPC0galxuICAgICAgICAgICAgY29tYm9zX2xpc3QkdHVuZV92aXJhbFtwXSA8LSBpXG4gICAgICAgICAgICBjb21ib3NfbGlzdCRWSUJSQU5UW3BdIDwtIGxcbiAgICAgICAgICAgIGNvbWJvc19saXN0JFZTW3BdIDwtIG1cbiAgICAgICAgICAgIGNvbWJvc19saXN0JFZTMltwXSA8LSBuXG4gICAgICAgICAgICBwIDwtIHArMVxuICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgfVxufVxuXG5jb21ib3NfbGlzdCA8LSBjb21ib3NfbGlzdFstMSxdXG5gYGAifQ== -->

```r
combos_list <- data.frame(toolcombo=rep(0, 64),
                          tune_not_viral=rep(0, 64),
                          DVF=rep(0, 64),
                          tune_viral=rep(0, 64),
                          VIBRANT=rep(0, 64),
                          VS=rep(0, 64),
                          VS2=rep(0, 64))
p <- 1

for (i in c(0,1)){
  for (j in c(0,1)){
    for (k in c(0,1)){
      for (l in c(0,1)){
        for (m in c(0,1)){
          for (n in c(0,1)){
            combos_list$toolcombo[p] <- paste(i,j,k,l,m,n)
            combos_list$toolcombo2[p] <- paste(if(i){"tv"}else{"0"},if(j){"DVF"}else{"0"},
                                               if(k){"tnv"}else{"0"},if(l){"VB"}else{"0"},
                                               if(m){"VS"}else{"0"},if(n){"VS2"}else{"0"})
            combos_list$tune_not_viral[p] <- k
            combos_list$DVF[p] <- j
            combos_list$tune_viral[p] <- i
            combos_list$VIBRANT[p] <- l
            combos_list$VS[p] <- m
            combos_list$VS2[p] <- n
            p <- p+1
          }
        }
      }
    }
  }
}

combos_list <- combos_list[-1,]
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


this function builds a list of all of the combinations that the user wants to 
test. 
In this case, we're comparing the performance of all unique combinations of the 
six tools.

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYnVpbGRfc2NvcmVfbGlzdCA8LSBmdW5jdGlvbihpbnB1dF9zZXFzLCBjb21ib3MpIHtcbiAgb3V0cHV0IDwtIGRhdGEuZnJhbWUocHJlY2lzaW9uPXJlcCgwLCBucm93KGNvbWJvcykpLFxuICAgICAgICAgICAgICAgICAgICAgICByZWNhbGw9cmVwKDAsIG5yb3coY29tYm9zKSksXG4gICAgICAgICAgICAgICAgICAgICAgIEYxPXJlcCgwLCBucm93KGNvbWJvcykpLFxuICAgICAgICAgICAgICAgICAgICAgICBNQ0M9cmVwKDAsIG5yb3coY29tYm9zKSksXG4gICAgICAgICAgICAgICAgICAgICAgIHByb3BfdmlyYWw9cmVwKDAsIG5yb3coY29tYm9zKSkpXG4gIGZvciAoaSBpbiAxOm5yb3coY29tYm9zKSkge1xuICAgIGtlZXBfc2NvcmUgPC0gZ2V0dGluZ192aXJhbF9zZXRfMShpbnB1dF9zZXFzLCBpbmNsdWRlX3ZpYnJhbnQgPSBjb21ib3MkVklCUkFOVFtpXSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBjb21ib3MkVlNbaV0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IGNvbWJvcyRWUzJbaV0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gY29tYm9zJHR1bmVfdmlyYWxbaV0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IGNvbWJvcyR0dW5lX25vdF92aXJhbFtpXSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gY29tYm9zJERWRltpXSlcbiAgXG4gICAgb3V0cHV0W2ksMTo1XSA8LSBhc3Nlc3NfcGVyZm9ybWFuY2UoaW5wdXRfc2VxcyRzZXF0eXBlLCBrZWVwX3Njb3JlKVxuICAgIFxuICAgIG91dHB1dCR0b29sY29tYm9baV0gPC0gcGFzdGUoY29tYm9zJHR1bmVfdmlyYWxbaV0sY29tYm9zJERWRltpXSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGNvbWJvcyR0dW5lX25vdF92aXJhbFtpXSwgY29tYm9zJFZJQlJBTlRbaV0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb21ib3MkVlNbaV0sIGNvbWJvcyRWUzJbaV0pXG4gIH1cbiAgXG4gIG91dHB1dFtpcy5uYShvdXRwdXQpXSA8LSAwXG5cbiAgcmV0dXJuIChvdXRwdXQpXG59XG5gYGAifQ== -->

```r
build_score_list <- function(input_seqs, combos) {
  output <- data.frame(precision=rep(0, nrow(combos)),
                       recall=rep(0, nrow(combos)),
                       F1=rep(0, nrow(combos)),
                       MCC=rep(0, nrow(combos)),
                       prop_viral=rep(0, nrow(combos)))
  for (i in 1:nrow(combos)) {
    keep_score <- getting_viral_set_1(input_seqs, include_vibrant = combos$VIBRANT[i],
                                            include_virsorter = combos$VS[i],
                                            include_virsorter2 = combos$VS2[i],
                                            include_tuning_viral = combos$tune_viral[i],
                                            include_tuning_not_viral = combos$tune_not_viral[i],
                                            include_deepvirfinder = combos$DVF[i])
  
    output[i,1:5] <- assess_performance(input_seqs$seqtype, keep_score)
    
    output$toolcombo[i] <- paste(combos$tune_viral[i],combos$DVF[i],
                                 combos$tune_not_viral[i], combos$VIBRANT[i],
                                 combos$VS[i], combos$VS2[i])
  }
  
  output[is.na(output)] <- 0

  return (output)
}
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## Calculate the performance of each pipeline
quantifying variability between rule sets

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzIDwtIGRhdGEuZnJhbWUodGVzdGluZ19zZXRfaW5kZXg9cmVwKDAsIG5yb3coY29tYm9zX2xpc3QpKjgpLFxuICAgICAgICAgICAgICAgICAgICAgIHByZWNpc2lvbj1yZXAoMCwgbnJvdyhjb21ib3NfbGlzdCkqOCksXG4gICAgICAgICAgICAgICAgICAgICAgIHJlY2FsbD1yZXAoMCwgbnJvdyhjb21ib3NfbGlzdCkqOCksXG4gICAgICAgICAgICAgICAgICAgICAgIEYxPXJlcCgwLCBucm93KGNvbWJvc19saXN0KSo4KSxcbiAgICAgICAgICAgICAgICAgICAgICAgTUNDPXJlcCgwLCBucm93KGNvbWJvc19saXN0KSo4KSwgXG4gICAgICAgICAgICAgICAgICAgICAgcHJvcF92aXJhbD1yZXAoMCwgbnJvdyhjb21ib3NfbGlzdCkqOCkpXG5cbmFjY3VyYWN5X3Njb3JlcyA8LSBjYmluZCh0ZXN0aW5nX3NldF9pbmRleD1yZXAoMSwgbnJvdyhjb21ib3NfbGlzdCkpLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgYnVpbGRfc2NvcmVfbGlzdCh2aXJ1c2VzW3ZpcnVzZXMkSW5kZXg9PTEsXSwgY29tYm9zX2xpc3QpKVxuZm9yIChpIGluIDI6OCkge1xuICBhY2N1cmFjeV9zY29yZXMgPC0gcmJpbmQoYWNjdXJhY3lfc2NvcmVzLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgY2JpbmQodGVzdGluZ19zZXRfaW5kZXg9cmVwKGksIG5yb3coY29tYm9zX2xpc3QpKSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJ1aWxkX3Njb3JlX2xpc3QodmlydXNlc1t2aXJ1c2VzJEluZGV4PT1pLF0sIGNvbWJvc19saXN0KSkpXG59XG5gYGAifQ== -->

```r
accuracy_scores <- data.frame(testing_set_index=rep(0, nrow(combos_list)*8),
                      precision=rep(0, nrow(combos_list)*8),
                       recall=rep(0, nrow(combos_list)*8),
                       F1=rep(0, nrow(combos_list)*8),
                       MCC=rep(0, nrow(combos_list)*8), 
                      prop_viral=rep(0, nrow(combos_list)*8))

accuracy_scores <- cbind(testing_set_index=rep(1, nrow(combos_list)),
                              build_score_list(viruses[viruses$Index==1,], combos_list))
for (i in 2:8) {
  accuracy_scores <- rbind(accuracy_scores,
                           cbind(testing_set_index=rep(i, nrow(combos_list)),
                              build_score_list(viruses[viruses$Index==i,], combos_list)))
}
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucGFsIDwtIGdndGhlbWVzOjp0YWJsZWF1X2NvbG9yX3BhbChwYWxldHRlPVwiVGFibGVhdSAxMFwiLCB0eXBlPVwicmVndWxhclwiKVxuYGBgIn0= -->

```r
pal <- ggthemes::tableau_color_pal(palette="Tableau 10", type="regular")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTUNDX3NkX3F1YW50IDwtIGNiaW5kKGFjY3VyYWN5X3Njb3JlcyR0b29sY29tYm9bYWNjdXJhY3lfc2NvcmVzJHRlc3Rpbmdfc2V0X2luZGV4PT0xXSwgYXMuZGF0YS5mcmFtZShtYXRyaXgoZGF0YT0wLCBucm93PTYzLCBuY29sPTcpKSlcblxuZm9yIChpIGluICgyOjgpKSB7XG4gIHRlbXAgPC0gYWNjdXJhY3lfc2NvcmVzICU+JVxuICAgIGZpbHRlcih0ZXN0aW5nX3NldF9pbmRleDw9aSkgJT4lXG4gICAgc2VsZWN0KE1DQywgdG9vbGNvbWJvKSAlPiVcbiAgICBncm91cF9ieSh0b29sY29tYm8pICU+JVxuICAgIHN1bW1hcmlzZShzZCA9IHNkKE1DQykpXG4gIFxuICBNQ0Nfc2RfcXVhbnRbLGldIDwtIHRlbXAkc2Rcbn1cblxuY29sbmFtZXMoTUNDX3NkX3F1YW50KSA8LSBjKFwidGVzdGluZ19zZXRfaW5kZXhcIiwgXCIxLTJcIiwgXCIxLTNcIiwgXCIxLTRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcIjEtNVwiLCBcIjEtNlwiLCBcIjEtN1wiLCBcIjEtOFwiKVxuXG5NQ0Nfc2RfcXVhbnQgPC0gTUNDX3NkX3F1YW50ICU+JSBwaXZvdF9sb25nZXIoY29scz1jb250YWlucyhcIi1cIikpXG5cbmdncGxvdChNQ0Nfc2RfcXVhbnQsIGFlcyh4PXRlc3Rpbmdfc2V0X2luZGV4LCB5PXZhbHVlLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb2xvcj1uYW1lLCBmaWxsPW5hbWUpKSArXG4gIGdlb21fcG9pbnQoYWxwaGE9MC41KVxuYGBgIn0= -->

```r
MCC_sd_quant <- cbind(accuracy_scores$toolcombo[accuracy_scores$testing_set_index==1], as.data.frame(matrix(data=0, nrow=63, ncol=7)))

for (i in (2:8)) {
  temp <- accuracy_scores %>%
    filter(testing_set_index<=i) %>%
    select(MCC, toolcombo) %>%
    group_by(toolcombo) %>%
    summarise(sd = sd(MCC))
  
  MCC_sd_quant[,i] <- temp$sd
}

colnames(MCC_sd_quant) <- c("testing_set_index", "1-2", "1-3", "1-4",
                            "1-5", "1-6", "1-7", "1-8")

MCC_sd_quant <- MCC_sd_quant %>% pivot_longer(cols=contains("-"))

ggplot(MCC_sd_quant, aes(x=testing_set_index, y=value, 
                                  color=name, fill=name)) +
  geom_point(alpha=0.5)
```

<!-- rnb-source-end -->

<!-- rnb-plot-begin eyJoZWlnaHQiOjQzMi42MzI5LCJ3aWR0aCI6NzAwLCJzaXplX2JlaGF2aW9yIjowLCJjb25kaXRpb25zIjpbXX0= -->

<img src="data:image/png;base64,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tVWrVjnLvvjii3LjjTfKbbfdNqYs7wv3p/yLFi3yTuY5AggggAACCCCAAAIIIIAAAggggAACCCDgS6BsArxve9vbHID77rtvAoR2tb///vud6To+7lRJx8xds2aNM5bpo48+OiH7Aw884IxZqXnc4RR0PM0f/vCH8p3vfEd27do1YRkNCG/ZssWZfsIJJ0yYzwQEEEAAAQQQQAABBBBAAAEEEEAAAQQQQGC6AmUT4H3zm98sy5Ytc25idtddd41xuOmmm5wbwSxdulTWr18/Zt5DDz0k9957r+zcuXPM9Pe///3O6xtuuGHMuL7aC/eWW25x5l100UWjy2jgWG/CpGO2/OhHPxpz0xm9mcO1117r3NDmrLPOcoLHowvyBAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQ8ClQNjdZ00G/L730UvnSl74k11xzjTzyyCPO8AmbN292nuudXT/72c86d4z3Wm3cuFH0xmW67LHHHjs6a8OGDaLDOWzdulU+/OEPy7nnnusEbbWHsI6hq4Ha8847bzS/3g3yqquukk9/+tOieZ588kl5xzveIXq3yQcffFBeeeUVp/wrrrhidBmeIIAAAggggAACCCCAAAIIIIAAAggggAACQQTKJsCrCOecc45861vfcgK8OoyCPjRpz14NrJ500knO63z+aGBWh1vQ8u655x7RXsCadPr73vc+ueyyy5wxeL1lnXbaafK9731PNGisgeGf/exnzuyamho5//zzneCvO6SDdzmeI4AAAggggAACCCCAAAIIIIAAAggggAACfgTKKsCrAKeccopzozPtZas3SdMbps2fP39CMNbFuvXWW92nE/5XVVXJ5z//ebnyyitl+/btzvALixcvlrq6ugl53Qk6vu71118veiM2XX9DQ4PoMnpDNhICCCCAAAIIIIAAAggggAACCCCAAAIIIBCmQNkFeF2ctrY20UcYKZFIyOrVq6dVVFNTk+iDhAACCCCAAAIIIIAAAggggAACCCCAAAIIFEqAbqWFkqVcBBBAAAEEEEAAAQQQQAABBBBAAAEEEECgwAIEeAsMTPEIIIAAAggggAACCCCAAAIIIIAAAggggEChBAjwFkqWchFAAAEEEEAAAQQQQAABBBBAAAEEEEAAgQILEOAtMDDFI4AAAggggAACCCCAAAIIIIAAAggggAAChRIgwFsoWcpF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AAEEEEAAAQQQQCBSAgR4I9UcVAYBBBBAAAEEEEAAAQQQQAABBBBAAAEEEMhfgABv/lbkRAABBBBAAAEEEEAAAQQQQAABBBBAAAEEIiVAgDdSzUFlEEAAAQQQQAABBBBAAAEEEEAAAQQQQACB/AUI8OZvRU4EEEAAAQQQQAABBBBAAAEEEEAAAQQQQCBSAgR4I9UcVAYBBBBAAAEEEEAAAQQQQAABBBBAAAEEEMhfgABv/lbkRAABBBBAAAEEEEAAAQQQQAABBBBAAAEEIiVAgDdSzUFlEEAAAQQQQAABBBBAAAEEEEAAAQQQQACB/AUI8OZvRU4EEEAAAQQQQAABBBBAAAEEEEAAAQQQQCBSAgR4I9UcVAYBBBBAAAEEEEAAAQQQQAABBBBAAAEEEMhfgABv/lbkRAABBBBAAAEEEEAAAQQQQAABBBBAAAEEIiVAgDdSzUFlEEAAAQQQQAABBBBAAAEEEEAAAQQQQACB/AUI8OZvRU4EEEAAAQQQQAABBBBAAAEEEEAAAQQQQCBSAgR4I9UcVAYBBBBAAAEEEEAAAQQQQAABBBBAAAEEEMhfgABv/lbkRAABBBBAAAEEEEAAAQQQQAABBBBAAAEEIiVAgDdSzUFlEEAAAQQQQAABBBBAAAEEEEAAAQQQQACB/AUI8OZvRU4EEEAAAQQQQAABBBBAAAEEEEAAAQQQQCBSAgR4I9UcVAYBBBBAAAEEEEAAAQQQQAABBBBAAAEEEMhfgABv/lbkRAABBBBAAAEEEEAAAQQQQAABBBBAAAEEIiVAgDdSzUFlEEAAAQQQQAABBBBAAAEEEEAAAQQQQACB/AUI8OZvRU4EEEAAAQQQQAABBBBAAAEEEEAAAQQQQCBSAgR4I9UcVAYBBBBAAAEEEEAAAQQQQAABBBBAAAEEEMhfgABv/lbkRAABBBBAAAEEEEAAAQQQQAABBBBAAAEEIiVAgDdSzUFlEEAAAQQQQAABBBBAAAEEEEAAAQQQQACB/AUI8OZvRU4EEEAAAQQQQAABBBBAAAEEEEAAAQQQQCBSAgR4I9UcVAYBBBBAAAEEEEAAAQQQQAABBBBAAAEEEMhfgABv/lbkRAABBBBAAAEEEEAAAQQQQAABBBBAAAEEIiVAgDdSzUFlEEAAAQQQQAABBBBAAAEEEEAAAQQQQACB/AUI8OZvRU4EEEAAAQQQQAAB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Gf1dFREQARFoEwE5eNsEVmZFQAREQAREQAREQAREQAREQATuEQJMghU+9y7smWk4gwMwybMsyjL4xTjqGJHpMJmW7wS+8D7c0YNcJ9HWxlmhCDJHfgyJvodQzZ5D1MozKLQGz0rBjY5y+cN0AOuVvq0HYQfGvUwvqk99ADEOCkTAaGvjnOeAQNUkWFMRAREQgTYT6LjeoMzslX/8x3+M7373u1hcXMTJkydx5swZfOITn4Bd76B3APXdd9/FF77wBVy5coUJL1N4+OGH8eKLL+LYsWObWgmyzabG9IUIiIAIiIAIiIAIiIAIiIAIiEBbCYTmZmFRK9WNRH3n7kY7M4m0jKPXWlryHb3O+ObvhBttH3RZLH0U6f770N/f75vI5/PIcuq/yh4lYAYCjESDKfRJQA7eVRb6VwREoK0EOsrBu8SO9nOf+xwmqI1kSl9fH7785S/7P9/85jfx7/7dv0OUnfJ2i3EU/87v/I6/eldXF7XsK3j11VfxR3/0R/jN3/xNPP744+tMBdlmnREtEAEREAEREAEREAEREAEREAER2DUCRjfVYkK1tURmm+zZZWSvvbSIENdXXOYmkLRYBERABERg1wl0lIP3137t13zn7tNPP41f/uVfRg91k65fv45f+qVfwksvvYTf/d3fxc///M9vC/LZs2f99Y1D2DiGn3/+eWoe1fDFL35xzc4f/MEfYGRkZM1ekG3WNtYHEdghgWq5iq99+XU4F7vRVUjBciwUE2WURubx4U8dR4ZThFREQAREQAREQAREQAREQATuTMBiQiyTFMvopm5VLPO961EKlxqrKiIgAiIgAiKwRwgwJWhnlLfffhvf/va3mUg0gV//9V/3nbumZQcPHsR/+k//yZdn+NKXvoQVjrRup/yP//E/2L97+Pt//+/jwx/+MPt5CxFq6XzmM5/Bpz/9aVSZUdU4extLkG0at9dnEdgugenpWXz196/gwDeP4szFQTx0PYUHp5J49Eov7nv9GF77/Rxe/d5b2zWn9URABERABERABERABERgXxMwmrugPINFyb8tC7/3pRqSXF9FBERABERABPYIgY5x8H7961/3kb7wwguIx+O34DVSDR/4wAd8iQXj5L1TKRQKvrPYrPcjP/Ij61avL/vTP/1TP6rXrBBkm3WGtUAEtkHARJK/+vllnLzWj6F8GLOpGt7tL+Ld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O/+YLJG/93u/59v7v//3/8JEAZtiln/605/Gz/zMz/gavI22gmzTuL0+i4AIiIAIiIAIiIAIiIAIiIAIiIAIiIAIiIAIiMBOCHSUg9c0/LHHHsMXvvAFP8rWJEkzyc9GRkbWOWPrkP7oj/6o/nHd71gshn/zb/4Nfv7nfx4XLlyA53kYGxtDKpVat259QZBt6tvqtwiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAjshEDHOXjrje/v74f5aUUJh8M4derUjkwF2WZHO9DKIiACIiACIiACIiACIiACIiACIiACIiACIiAC+55AxyRZ2/dHUgBEQAREQAREQAREQAREQAREQAREQAREQAREQAT2HQE5ePfdIVeDRUAEREAEREAEREAEREAEREAEREAEREAEREAEOoWAHLydciTVDhEQAREQAREQAREQAREQAREQAREQAREQAREQgX1HQA7efXfI1WAREAEREAEREAEREAEREAEREAEREAEREAEREIFOISAHb6ccSbVDBERABERABERABERABERABERABERABERABERg3xGQg3ffHXI1WAREQAREQAREQAREQAREQAREQAREQAREQAREoFMIyMHbKUdS7RABERABERABERABERABERABERABERABERABEdh3BOTg3XeHXA0WAREQAREQAREQAREQAREQAREQAREQAREQARHoFAJy8HbKkVQ7REAEREAEREAEREAEREAEREAEREAEREAEREAE9h0BOXj33SFXg0VABERABERABERABERABERABERABERABERABDqFgBy8nXIk1Q4REAEREAEREAEREAEREAEREAEREAEREAEREIF9R0AO3n13yNVgERABERABERABERABERABERABERABERABERCBTiEgB2+nHEm1QwREQAREQAREQAREQAREQAREQAREQAREQAREYN8RkIN33x1yNVgEREAEREAEREAEREAEREAEREAEREAEREAERKBTCMjB2ylHUu0QAREQAREQAREQAREQAREQAREQAREQAREQARHYdwTk4N13h1wNFgEREAEREAEREAEREAEREAEREAEREAEREAER6BQCcvB2ypFUO0RABERABERABERABERABERABERABERABERABPYdATl4990hV4NFQAREQAREQAREQAREQAREQAREQAREQAREQAQ6hYAcvJ1yJNUOERABERABERABERABERABERABERABERABERCBfUdADt59d8jVYBEQAREQAREQAREQAREQAREQAREQAREQAREQgU4hIAdvpxxJtUMEREAEREAEREAEREAEREAEREAEREAEREAERGDfEZCDd98dcjVYBERABERABERABERABERABERABERABERABESgUwjIwdspR1LtEAEREAEREAEREAEREAEREAEREAEREAEREAER2HcEwvuuxfd4g6PR6B1bYFnW2jrbWX9t5U0+GHtB7YTDPzjFbNsObMdsWy/GZtD6hEI/GNOIRCJo/Ltufzu/W8G40UYnM248dtthu9E6rTje5lgHtdOu89jzvI2ae8dljeetOY+D2mncUVA2jcd3LzJubONOPjden4Zxs6VV1/heYHz7+dfIaiecGrcLev417m8vsGm8HvZCX9XIJyjjxuPUyYwbj10jt518Dsq48Zrai4yD9jGN504r+ipjLyjjxuO7Fxnv5DxrXLeRcVA2t9sLaqfxWWkvMm5k1djmO31u3C4om8Z97EU2QdvVyKaZ96o6H2OvFXXpZMaN97I6t538boax4Vove5FxK/oqc/4FtVNno98icFcJ8ARW6TACP/7jP+49+OCD/o/ruoFaVygU1mx8+tOfDmTDbPTtb397zc6v/MqvBLbz3/7bf1uz87/+1/8KbOfnfu7n1uycP38+sJ2Pf/zjvp3HHnsssI35+fm1uvyDf/APAtv5+te/vmbnP/yH/xDYzn/+z/95zc6XvvSlwHb+xb/4F2t2JiYmAtv50Ic+5Nt57rnnAtsw+69fCz/zMz8T2M6Xv/zlNTuGU9DyH//jf1yzY45b0PIP/+E/XLNjzqOg5fHHH/ft/NAP/VBQE565juqMf/Znfzawnf/9v//3mh1zvQctv/qrv7pm55VXXglqxvvMZz6zZiefzwe2U2fzqU99KrCN119/fa0uv/iLvxjYzuc///k1O3/wB38Q2I6pQ71db7zxRmA7hkndTlAj5tjUbXz2s58NasYz50rdjjmHgpb/+l//65odc04HLeZaqtfn/fffD2rG+9jHPubbMdd60DI3N7dWF3PvCVoa+ypzLwxafu/3fm+tPubeHLSYPqHO+Nq1a0HNeB/84Ad9O6bPClquXr26Vpd/+S//ZVAznum76236/d///cB2/v2///drdprpq8yzTb0+zfRV5lnL2DHPXkHLe++9t1YX8ywYtJhn0HqbmumrzDNx3Y55Vg5azDN63U6xWAxkxrwr1G2Yd4igpbGv+qVf+qWgZrz/+T//51p9/vAP/zCwnV/4hV9Ys3P27NnAdn7sx35szU5QI7lcbs3GT/3UTwU147388strdn7t134tsJ3/8l/+y5qdL37xi4Ht/Kt/9a/W7DTTV7344ou+nSeffDJwXWZnZ9fq8o/+0T8KbOdrX/vamp1m+qrf/d3fXbPz53/+54Hr88//+T9fs3P9+vXAdp599lnfzvPPPx/YRmNf9bnPfS6wnT/7sz9ba1MzfdVv/dZvrdl56aWXAtfnp3/6p9fsLCwsBLajDUVgLxD4QXjlXXUza+etJFCtVmF+TOFJhsbR1e3ux2xXt1Gr1ba72br1+NC4ZsdxnHXfb3eB2bZen2bsmLbU7Zi6BS3GhvkJwra+z0bG9TrVv9vJ70Y7zRyrdjA2dQta6ow7lU0rzj/D9m4zbrzG98L513geN8O48V4RlHHjtbkXzuPGY9Wq+2hQNubcrV/j9fM4yP20kXEz51+r2LTj/GvmPK4zboy6Mbx3UtrBuJlj1ch4L53HzVzj5njUt6//3skxqq/bKjaNdlp1jdfrGOR3pVK55X4RxEanXuN7qa9qZNzMNd5op5lrvPE8bsV9NMh5V9+mHffRvcamFfeKxijwOrvt/u5Uxq24xg3D+vNAM31MI+Nm7Oy1a7xVjLd7rmo9EWgngR/E2bdzL7ItAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiLQcgJy8LYcqQyKgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIwO4QkIN3dzhrLyIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiLQcgJy8LYcqQyKgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIwO4QkIN3dzhrLyIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAiLQcgIWsyEGT3Pf8urIYCsIfOMb30Aul/NNfeITnwhk0mRm/cpXvuJv293djQ9+8IOB7CwsLODb3/62v+3hw4fxwAMPBLJz6dIlnDt3zt/2oYcewqFDhwLZ+f73v4+pqSl/2w996EPo6uoKZOev/uqvUCwWYds2Pv7xjweyYbJCf+1rX/O37evrwwc+8IFAdmZnZ/G9733P3/bo0aM4ffp0IDvvv/8+zI8pZ86cwcjISCA7pi6mTqa88MILSCQSgewYNoZRNBrFiy++GMiGOUbmWJkyMDCAJ598MpAdc86Yc8eUEydO+D9BDJlz2JzLpjz++OMYGhoKYsa/psy1ZYphYxgFKeYaN9e6OUbmWAUp5l5j7jmmmHPGnDtByvXr13H27Fl/01OnTmF8fDyIGbz99tu4evWqv625psy1FaR885vfRDab9Tc117i51oOUL3/5y/5m5l5j7jlBytLSEr71rW/5m46OjuKRRx4JYsbnYviYYu7F5p4cpLzxxhu4ceOGv+mzzz6Lnp6eIGb886aVfZWph6lPkNKqvurixYt47733/Co8/PDDOHjwYJDq+PebTuurZmZm8Oqrr/o8zPVtrvMgpVV91Xe/+13Mzc35VWimr/qLv/gLPzt5M31VoVDASy+95NdlcHAQTzzxRBA0/vNNK/qqd999F5cvX/brYOpi6hSkmOe/Tuurrl27hjfffNPH0Uxf9dZbb2FiYsK304l9lbn3mXtgkGL68Hpf9eCDD2JsbCyIGbSqr/rrv/5r5PN5vw5B36tqtRq++tWv+jaa6avm5+fxne98x7fTzHtVq/qq1157DdPT03599tJ7VX9/P5566im/Xjv9p1V91fnz53HhwgV/94899hiGh4d3WhV//Vb3VbFYDB/96EcD1aVVfdXk5CRef/11vw4nT57E8ePHA9WnVX3VK6+8gsXFRb8OzbxXBWqENhKBFhOQg7fFQGVOBERABERABERABERABERABERABERABERABERABHaLgCQadou09iMCIiACIiACIiACIiACIiACIiACIiACIiACIiACLSYgB2+LgcqcCIiACIiACIiACIiACIiACIiACIiACIiACIiACOwWATl4d4u09iMCIiACIiACIiACIiACIiACIiACIiACIiACIiACLSYgB2+LgcqcCIiACIiACIiACIiACIiACIiACIiACIiACIiACOwWgfBu7Wiv76dcLuOP//iPYTJVmiyKJqOjyQZvsqXuJHu6yeb4hS98AefOncPKygqs/7+98wCzpqgSdsuaRRRERUQEFxFEMS/msBgxbnBdVgyomPOqGAFzXF0zKovogyIIBkQQRTAQTICAigqIKEkRFDH76/z37c/3cuZ81ffe+WZGBr5TzzNT3VV1Tp1661RVd92+fa90pe4GN7hBxy/k8quMN7vZzWZCgZ6Pfexj/S9M/uEPf+j49VV+9ZJfab797W/f8Wucs9jW0nPlK195/IvP66+/fsevR9/3vvcdtE02/Boyv6r8l7/8pcMmufz1r3/t1llnnd4+fpn9Hve4R9M29Xz5y1/ufzl4bm6u4y/Lw+oLX/hCn86vx++xxx7d1ltvPeYWGfMr99gCZ/hQP/IHH3xwx6+zk3/1q1+9e9jDHrYa/xYbGfNrntRNe69xjWv08nA77LDD5tljm2SDkaShBz7oOeuss/pftkaPNmZGWQ/tgSl9TcyvYt///vfvDj/88O5a17pWd9Ob3rS74IILuo022qjbbbfdxmzUI2Nk+bXeDTfcsP+VW34ZFpZ//vOf57Ut+2ZkrB+jgzaceuqpvQ9c85rX7K573ev2/oNfZjYYFRljGwz51V+O+QVfAjrxg6tc5SrNsWKbZCwb2PJLy7SJ4zvf+c69/s9+9rM9I8baNDbZfn5VFt9hvNF//CJwnge0JzJWD7w/8YlPzJPHXvrvla98ZdOP8Q/6lF+bhg/2O0bxQdhc+9rXbrKBn32FHlgSmH+03/qj/wyxcf5jzmIetI+o3/6PuofYRD3YxS8s4y9wxR7HgWP9Rje60by+im2CDfn8oja6kIWJjPC9oXks6oEN/XPuuef2PkObnE/jGGPu2WNgzhlijD3HHHNM30+MS365mXk1ztH6TWQD4x//bY7ZfPPN+3kCf0Yef8aWzKalB7uQo1/hTB8xT9NG2tZaMyKb7H9xjrn+9a/fzzXomzQfYwNyW221VcdYZV7HJgK+rY2O1aE5J+rhl8v/9Kc/9fXGccD6x6+Ib7bZZvP8JrNhHWe+Y0zTL6wJ+g+x/Z/X05aeVh/bR7R1FjZxjMuYOQzb+GVq+g+dmU0P8W//mB/22Wef8TqAPHN0XOsYG635GBW07aMf/Wh3wAEH9BqZ4zIj1kzHaF6r/mbGWM9+++3Xryd5PXAeHlqr1IMfooPxA2PmLecI1wrHSGsdj3re+ta39r5vH7fmiKG1Cj32+6GHHtqdffbZ/XzqWoBN2EYYWqv6zNE/GcsGewjE/Io68x62ySjPx33h0b/MJq4DlMlrTZ6PJ+lpzeeuFXmtinpkrI/YR3GtmHTNK+PPfe5zPeM41zmO7ne/+42vdYbYyFg/jnNMnAddB6exOe644/prxbve9a7jtcY2el3Zmo8jm//93//tfvjDH3bM5euuu+54rWHtwY8mrVVDepzPbYd9xHreWqvQM4kx8zi+jX1Da5W2ZMbMFZPmwbxWqYeY+e0FL3hBvy4xV3qtxbUo/T5prYp6Pv3pT3fve9/75tnvPIzf0aahtaqlJ9ffuq4cmo9p03//93/3bTvvvPP6eTCP9aG1KtqCnmc84xn9PO41n3MNvNE5tFa19LDusta1/Ncxlq9xop4hxqyjv/jFL2Zeq9TjHMoYjdc8zOeT1ipsiozj9QRrwvnnn9/P0fCZdO+pHvwvX8/EuWbaWhX10KY4jma5r0LewNz/4Q9/eDwOuM7Bb/M6OrSORz1cD7jW6T/T7quUJ3ZswqZ1PTLLWhX1OP5gyxzlWNcPWvdV0R6PM6Mo7zw4tFapw7g1xuJ6OjTGlc+xtkX/GVqrsmw8b82xQ2tVlJt0bH/m/YFJMpV3xSNwpdHiMXfFa9bCWsSE+PSnP73fdERygw026C666KJeyT3vec9u99137xfXaVrZIH77298+WIwJ9w1veEN/sz1YaJQxTY+y02ybVQ/6uBB44xvfuJptmY11T4uzbWuqh3re+c539ptsHC+kTZQnsBDg5pH/muhZpe1SexbTJnTJiMk4+p/15JiNGi7OCGzaXnzxxf1G/zve8Y4+bTH2LDUbDFoMY+1hwZyFDRfDXFTEwIX1YthwUcMFBWGhfYVMlOecG2NuZgmRDTdaXNzxN0uQDZt2UQ8bRr/97W/HKqK/kBjPJ7HhAhOfnBSiriE2s+ixjtvc5jbdu971rv40tgk2v//971frW+VyHOexqAc2tGnacme7huacaYyzPUvNJo/xhTDGNu3hRs21asj/nGOMnUeH2AzpkalsHKtDfjykR/kYb7HFFt0HP/jBPimziet4lBk6lk2ejxeiZ6Fs8Fc+oDQ4Z0Q25hGfcsop3bOf/ezxvJTl7SNloj2kwehpT3taf2PKufVxPCnIhk0GQtYzSVb/iXOO5eMYNW1a3GKzJnpabGZZa7J92Z6FsFHXYtnEfsx9NQubPEbjWqWNs+ixrHGLzSyMnSPQ02KTGWf7rZ845q0Jm6iL47hWmbcmbOJaNYue2A7Ke75UfpzZZMbRx7Q3xy02lGHdfclLXtJ/iJPnqKzD82wP6WzCv+hFL7LIavF6663Xf9BlHZmNAtP0WE4/zH5MPm1iw5APM2cJca2K5dHzuMc9rjvzzDNj8sTjFpuF6lksG/2hxQbjM2P9dahhLXtajKfpGWKj/7ke6SP2MXaZ15pzyMce9XA+LQyxyev4ND0tNsgsVE/LntimfD1hH0f7ZmHTYquOqLPVV5YbaluUp2xrrVKH8Sxjo8VG+RxH26L/DLHJ8p5PmmMnsVG+Fcf+XKg9LX2VdvklUK9oGPXdq1/96n5zd7vttusOOeSQjk8deXqWT3K/8pWvjDeJJnUzA57NJJ7sIBDzydKznvWs/pxJgI0KLgj41HEoqIcJksAFC4FPFmMgfZJtQ3qiDo75BIvATWbLNtmw8BloSyvE9GybepBjgm4FFu5JwTbJWEY8QUIYkueik5ti+fPkJX2l/BBjbSGfyd/g5r9tko1xLKsMMZt72GGQUdaDH37gAx/oNt54Y4v2sZu7nLC5m4N6SIfRa1/72v7JuVwunl/vetdrspGxfrzJJptEsf4YhjyxEYNs7CsZ88TS6173ulh0fGy/Gcex8vKXv7wfm7KFza677jqW9YCnBoe4U6bFhqersv3qI+ZpgzwPLEQPT/HoW+g76KCDiPqLMtixWUIf8bQYH/wM2c9Tiq15RD9Wzy1ucYtev/9Y6NksM0T/MY3YNsH2TW9602qbu8wRuf/RzUVnnCOn6aGuobHKkxAE/cY28aFT3rjvC6Z/eR7LbHi6A5tjePKTnzxxjFE229NirN9G3XwoOI2NY0O5qAd7DbLhfBbGyhHjV9hBcM3gyZ1Z/M85xjjzy2zw7+23376vK/7LfoefteZj+xw9fHiRA+PgqU996rzkM844Y7yeRjas44zdlh4UTFpPsx54xbVNA3jKYihkNkNjPG7uouvNb37zPDbxWuHEE0/sXvrSl443dymf5enrOI84H1OWQFuiP22zzTZNRq95zWuaa9UqLavrYa4bmk/1H2WNYcSHDHkcOEdYzjivVbJRj+WGYj4MiSGzsd+jPXFMRlnT41qlPZkxclzz8OTurME2RVviOpD1tNYqygzpac3nca1wrbIe9XhOvMsuu8TT/ph5uLVWyUbG8RqwNY6mzfmZcZ5jomGttYp822RfRhmPeTqTD0RiiOOHdPXEMvmYOuz/vFbJZpoe2hH7aKjNszLWxrxWmZ4ZM8anzYOZDbq47mYd5wl9Ql5HSGNt4po3r1WsUwbnP8+No//wLQZCqw7Lt/RQf+u6Mq9V9pVtam3uDo31uFZpC3pe/OIXNzd3fRLYssau45HNkB7m8ey/6mnFLTaU45swOcyyVkUZ/DXOZ7Rj0lo1xDj7/dB9lXWrR/9zPdJH4lxjnrIxznpiXj4eWqsoB+OXvexl89Zx5YfW0bxWqSdfD6iH2HnNuLVW5Tbl64mnPOUp8/oo6o/HWU+LreWH1irziacxwncMea0y3RjbWmMsrqctNsrnONsW/SeXnXae59hZ2EzSmfthUtnKu+ITmLyjdsVvf/9VLBZpFhoWQzeSbnzjG/efDDHw+YoEX4GcFD70oQ/1FxV8RYjwmMc8pv+a5X/8x390//7v/95vUrA5wNeBPvWpTw2qUg8TJAs8FyzYttdee/V6EGSRZCJm4h6yraXHid7NA/Tw9QDS+cu28ZVo2GAHC5/yPPYfg/JMdG6AoFvb1KO8F7ro4BNDGJOXF+5YB8e2ScYwgs3//d//9a+FGJLHlkc+8pE9P9r4nve8p++rFmO+6p0DfcDEaTj66KPHfiMbY+zZc889LTov5mv6fGUbfyBgF1/ziIyRZ9Fmg4Kvks8aMmP8j5sBbhpagboJXDhoU2QjY/TwlarWhTsbfHwtLQbYEOwr+4iNTL6O1Ao77rhjz4T+w3/0I+w56aSTxv4nG5/Yi7q++93vzuujmNdic8c73rHfAMr2RzmO4zxgX+nHsJmk57TTThs/TYIu+pl5RDY77bRTL/8///M/3fOf//xB+/HJ1jyiH6OHryWecMIJVDMep7Dna3qTgmxgy/y39957r1Ycm3P/o/trX/vaeI6UjX30lre8ZTU9JPAVVr7ynccqfthiw4dkQ2HSPBbZMA74On8O2DxtjMW+GmKc20I9fP2ZJwuY2yKbyBiGMUQ99Esco7CJfcUcMcQ46mT8a4drBmOUumf1v6gvHkc2bJ7xFBJz/lCwPfRFa87BHvXgWzlw45zXTtpBWmQDY8Ye49MxEXWxkTy0nsa+gjGbDmxGty7imc+HQmTDHDFtjKuHNsa1irbxRDF9jR088TEpnH766fPmEedjZLgROv7443tx5zC4tRhxHaQdCNB3rudRT69s9I8bLjbUp82nlid2gyKOA+xiPW+tw3mtkg03STnoa84R5OdvJUQ2MGB+JkR7eHo7B69V8lqFPS02yD/2sY9tbmZn3Z632NAnfH12KMS1yr5q6fGrxVEPbWatkJtrFWXg1mLM66JyQMeDHvSgfj2P6zhsImNssK5J4yjr5zwy1o9b5dRP2+Jaldkw7w7pOeqoo7pPfvKT89S7VpGY2VjnPIHRyWajD6fe//739+tB65p7kh5ts4+y7ni+EMbaGtcq2bQYM6cOzYPqimywi9cRMQ8zr04KbGzxIX1eq7CH1zs5//mNqqhrVv+B8ZAePvBrXVfmtQo/pk2M56E2kcfrFnKg/5A3qOfYY481aRzDEz/h+oEQ76tcx+0rGbf0sCHNa1ymhUlskKXvc5i2VrX6Cgb6Suu+yjpkM8TYcsSt+6rIZlJfRT2TjmU8zR7bNrRW6X/cU7XC0Doa16rYV5OuBx7/+McP3lfN4sfY9+53v3ve9UTL5lnZRNnWWkV+bNsQo6H7qqjfY/2oNTaG7qviGFUP8Sy2xfLTjltzLDJDbKbps634qGvGNJnKv2ITWOs3eL/0pS/1Pcz7EtlsjIGLe5485KvTLBhDgYHvzcE555zTF3vAAx4wLu6xeWzetRa/qAdhn+DUNvWw6GPTTW5yk6ZtQ3rcQOBGmhD1mBdtk412WMZ29EpG/7Aj56kbburJZZC/973v3b930Dx15ji2KdYPG97TyN+0IL+4qWPb0HPEEUf0m5nocaFWJ+87Y+IlsFgceeSR/bHyxuhhY6AV2LCEh3bASD+I8s985jP7p8izDVEnzGPIjKmDp18J1OOGLedM/s997nM57G9mo02ykTFfSWdTrhVYaNnQi5/uw4YLN8cDcrLlfYKtgH0ysV5jyvuULXr4eqfv4eKmKQc22g1cTBJabJ70pCdNZawe5wH7Sl/FZvVYNsb0H+8A5jUTBPqfbwfIZpp81MVxZmRfxT7iZs0PqbJ861w20/oI2ayXuYILUuZI2eQ+ynXy4YdjlX430FeZzROe8IRxX1Mulud80jwmG5g5DhxP+Cwh+lifkP7FOWdNGGe/yYx58ij6K5/eG2hrHqOxr+I4YG4aCrGPWDOYwxwXs/qf3OL6mNk873nPGz9J27KFp/riOGjNOTCepIfNS25Yeeo0BtrofAxj/As9bAy0Atxpi+MJ1q6n0Y+dh1s6TPObDp4TRzazMo5s0aFttO2JT3xiPzZIHwrcbBOY25lz4lpFm3gPnvM4cwSbcATmpElBO2SEDvVEOeZ2xm8r5KdFKUO7WusB7zKP63C8UcFvo9+gY+edd543R9gf+rhzhHbRdgNrlf2d5wjL8L7JGODq/O/8YcxGSosN8nB84QtfOB57USfH2svxEJtJjJV3zsGf2UxoMfbDujyfMy+y0UbwWoXjPA/7FGlrszn3kWxoE0+tEfClf/mXf5nX5j4j/ONDmFbIfhyf5MrzINeWBup3rcpsmA+iHmRcI9gwzXOObaRcZMMc51xJXgx88BHrj9fM2Y9ph/1JH+V+inrzsX48C2Ovq9AR/SaO8VmuJ3j3KE86EyIb+mrSPNwL/O2f4zza4f3NpGusqCMe5/mUvKwn9hXvrG9dVzoOnAcZ45PWKOphHYu//UCagb52Pp6kx+sCrq0IjiPnYdnQVy3GlCPw2quhdbAv8Ld/mU3Mm3YsG9qW1yrtUAdjjTnekO+rJrFxjVPWeSj3kYxk41PXyk2K830VZWf1Y17BkcdBi030u2yLczwfdOf7KteqWfuKuq1f/zGe5sfZRv0Re52bOJ6VjeMbGUMc66zrBNuW64/yHOf7KuXVTaxt0/o/M3KMRl0cD9k2tFZl+XiObV4rMMfygE8MLTYxPx/HtvKtg6wvl6/ztYPAWr/By5N/BAZFK7jAnnzyya3sPo0fnWLScxONRcINOwrwyTRfy+LpBiZunoZtXRyrx8nUT8u1TT1O9N7MZNuG9GALtvEEHfZkPSyg0TbZaAfytlEbSdMO5G2jurFNPS152pafCG5NmLYp1488n5bCO19MUF8M8Wv82mzb0MNFHoHFJS5ipPFpuheNvOeUr2kQlDfWnj6z8Q8euR+jHuT5VJZFB1+xvVlVvtiJjPU/NhoIbFTEp6a5QdVe8rXJmxoWFwJ64BD7Iy+8fB1whx126MvzDza8egI5fSSy9euJUQ8/qiQT/YfYRd2vDUW23Jw8/OEPH9fLAU85uLHBuXItNrMwRofBecBzGaOHevUn84nZ0OIF/DHwXjLYZHn6mfclxZD9OTKSLWMn9hEbf/Fpc/TJMer2WDatPuJrZjE85znPGd/4otO5IrKJfaR87Gu+pshYZSzhg7aDejIbxjVB/4vtoN1D85jlKENdjgNvaPVBdOPzcG0F55w8DlqM401L1JXZOMfw9ABP18UNNux0XMMmj9HYV7KhTVtuuWVzYwH/yX2kTtoNG/0X//vn0Y+AWn9sQ2uDLrNBF/KtcYBOnkCPIc85+jF6HAfOt8qhh6/avfe9713ND52P9T/0aEtuk091Op5cqyxPfehxjtCf4uYgZXgylzGRQ2QzjXEe4+rSNvqPvmaOafWRtrnpyGYLc47smI+51kCH4wj/tbxjMzPKdsiI8aseyjjGydfGeGNKGeanHBwHjkXrZxPMPGzDzjhH6DfIwebCCy/sVUebrItxyxwRQxynsiE/zxFRJq7BrpHYG9cq1mrm3ehDbsxhv36gLvvNelyrOLf9ssmMYx0teecc/JCQ9djn0XexB/vZ/DTAmuBcI2PaalC3OkmPfSQj+krGXI+wTtvnjKPYx+iwzRzHkP1YndiW50E37FrzoGych9VjO2wXdbNpFecc0jIb5mHmOMdVth+2jEP7xj6EQfTjPJ8/4hGPWG09Z64cCrajxdi2DclqWxzjrbUuz4O8AiA+cR/Z4PfMRWzk5T4emvu0Q0b4ENfDjhnrV57rVcsOtY101zrL6n/kOb/BP19XRn+mH/khqrzWaQttZL6j3fxwNUG7ielr52PK5LXOPnLOhBu+6DjSTtvg9ZTntsk68UP5uyb0RqV/Q2ud7aJ47OMo3lqrtCeOA2zjmgeGhnxfFdlkxl67yCjOQ7GPZCQb9NCv2qS8NsQ2ys88Yuac6MfKq8+yfMMmMoo20e/o8b6OtbIVsJ0y++yzz2r3VbAh6MesJfqPumyLa459o/8YRz+GjXK2X/9jveSaq3VfRZ2ZTR7j2pU/QDPdse6coR/St5FRlPcaRx3Eysc0bZOR/m8bLZsZOUbNN862ySz6uGWnxdjmesocG/1G2czG9FZsW+krXp3S0teSq7QrNoG1foPXT7TihX/sctN/+tOfxuR5x+pwwlcmFjLNG4SWPvW4kDGhEJTNxzwpSci6hvRE+ahTPS5c6lOPdiBvG7WRtCxvG8lDl3o4z/LY4Q9POfG2Fgl1tOR57+NHPvKR8Y2MFzfUl4OTvPbbNuxgoSNYR5SViWk+Iai8MXr4StBQUE/kT9koz1fFeBKDr/nmJ5CG9MqHfHV7k8KTrnGDlzLxySRtkg35BPTA1k/KSXNRs6+4MbU+8gk+NSRj8mXrQq2vUd761aP/eJER2ciWi3zLo4PAZlnUy8UjocVmoYxzXZ6jh/ecyrqv8G//vNmLafpNlt9///07fgk6htgW05WTP2Ml9hELu+9A9AKCG5KhIBv05j7iSZQY0G2aY4y+0ybKcmwf6U/2I/n0EWNVm7SRvMxGPXLUnyirXKw7z0POI/YNF7cEfZhjvkr5qEc9isPVQmQzxHg1oZQQ7eNYxmya0r8yQgw28f2jeYxGe2TDOMA2+yNWr4/EPuJVDwT7RP/F//CbWL8XruqJuqMt1E+fIi9rympT5K0Ox7zt14/R4zjwBs5xALMHP/jBPbfsh9F3tIen6wj6j3XH+Sz2j/5DOdKdI5yz8rrAKylaF9KRDbomMbZtlMtB27gBYI5p9ZE2eQPmuIi6YJ39d9ttt+2L6I+ZUZTXDtLw0cjPfuAmSRsf+tCHRvFeZl7C6MR6ZWv9+Jp5ti22Sb/RJn6IhKB8nCMsw825IfYxaerTb50jLE+sHRwrb7+5VhlHW/V7xw9+4HpuHjoJrlUc237ZZMbaSllDlLfdbhpnPdaN7+q/jlX61yAb5xoZK085dcc05bRDNo51rkfwR+caxhGbIzHIN6ZxnP343ve+d1+EPiLPdpDo/Cl/7NIm2TgPq6dXNvoX28NTzdhjP1DGNsrGTUfrj5tY6oz1Zz/Sj/N8Tl/n9dy5W70x1jdajO2HWD4eyyaOcfzD+uUYxwPyua9kQ38wn++7777dK17xitX6WH+KNnCsHTJi/HI9rK9Zvx+W8DSf33bJuuK587CM4lzhpnnrutL2aBeb53mtkwGMbPfd7na3vnrz9A30WSavdfatm5TUab2x/2RDX8E4r3XWCXf5x3k7cuF4aK2TMWXkxnEO2uhaZdk4jljjsUt/Qods1RfZZMbWYdmoWz2xjGzQ89GPfnTsI35g7XwU26juGNtXcrQfZGzZfJ5tgo33dayV0X51sI5ShjUgtoV89enHk+4ZHKvIqUe7jfVj2NBGgjZ57bzH6PVHXnP1BUb/4lqT2TiPK6/PG6vDWNti21qMonxcY9WjvOfE2pbHmG2MZbVDNi19jhHuy+k/f8ch93vUO3SMbY5H1+BcVptatuSytjX3VS5X52sXgVXf61u72jyvtTzFQXAwzcscnXhRYbmcz7l5TjwtXepx4lMm6jPNCx8vQqM+9SDnxbhy6vI86yFfXS09WV492kG+bVQ3adrBMcE2cowO9XCe5bFHXU6UXsBR3qCOlrwX19bLBK5O5Y2V90LKtmGH9aMnPsmCLPU7IXPuIqe8MXpaCwgyBNsR+ZMe5d2MJfbGnTKTgnopYx9HnV7MqMM8zpWVjYspemDrUyGUpW+4+LQMsvInn+DFp30Q2cY6zLd+mdiPq7Rd+nRFZMux5S2HnthHXsCrn3LIERbKONeV9URGfQWjf9bvObHcszx5+h/HBBmvOlv1Xzv0MXjGPkKvtnCRhx/LOerxWDbIWb99FG9oKU8ZL4C0DfnY/5SJtlmPMeMGu2xHHOuZjXq03zGLLv1ZPaQ5D1nOdqjXWL3IUH/UQZohshlibFnZeW4cdUfGHBO0mWPYtHiQhy3RHtuAHmyTEWUNlol95I2jjBwHyKAntiMeq9M42mL/R0aUw0ewS19Rllh5+0g/Js96vUnhHHtjG7MfOh/LAz2237aSRohzROyf2BfokU3Lb1Zpav+3bfaxeiidGbfYqFXb8InIWD6Us42yUjbG9r/zAnZ5rK7MKMprB2lwVpZz54jYf3nNiv2GDMF67X/rpx3muQ7HMSFbbdKfrUM91CF/bqD8ACD2MWXUZx8bk2fQRs6zvGuVcbxusE+M8YPMxjriWmH7rTczjvxb8rLR1qxHVvBhLPCuQ/3QPPTKRh8zz3PKqJs+si+U0w7ZaA/15nFgX6FzUohzDDLaDWPytBEdzgnai12Oo2gLZdWj/yhDnu1xzolpltN+5aMfUJ7Qqn9VzvxrnNgO+sfNMto4bT3Xf1uMbbN15tj+imMcPfqbbG1zlvdcXrJ2/pOR5WTlubF2aK/9aNus33HCedatrhhrh+2J9asLPdavrO0xHTtom3oop7x+SJrlOSboY9EPotyqUqv+yzra4ziiRGRDnTKxTeol3XZH/411cWyZIcaUcR7iOAfb6lolG+2iPHYSYl/BIq7Hkc0QY9sYdec+oh78JvaH9WoHZQiyWnW2+v/sx7KfJpdtch2TtX0ca3Q+JU2m5qtPedKzDZ7HNqpH/zHWj9Gjb8pUPZ5TxhDXqsxGxsrbxpYe9GlbbpttnCavTcp7Tqxt2iKTli3aIZuWPm0yVm+sc9bj6Nswc7xEeW1q2RLLcWxbc3qdr90E1uoneJk8XNC8iMru4OLjRUbO59zB6cTR0qUeJ3DrjfrU46TmRWLUpx7kXGiybUN6kFFXS48LJ7ZFNtqBvG3URtK0Q3nbSB66YluzPPZoi/JOxMgbbFNL3jLKWcb0GFvGumwbdliHm0dRLjO2zcobRz1R3mP12GbTo7xpxLlczPM49hVp2BDTOLe9ymg/59okGxeurIeyXqR4QUDfZhvlqI9EJvaN8rF+9eg/2uyNY9TDseXRQbAdq85WPdUVOZCOXA6TLlwt69OEnkc9uQ7L2H7PieUe5c2Xm+f2g+fEttk8WMX6o177U45RD8dZzvrto3xDgG7T7D/6P7KhTNYT67WPbIc2UiazUY8cYzv0Q/VEectlNo4x24dMy39JJ1g/bcqsVpW49L/9cWnKqqMhNuq0zZSWjTpiHnbyR4iM1WOblSWWbWyjY6lVPspyPKlMZKNcZGQacYuNbdVGxzzl1a2PKa8fUCb7oazgYVBPbkfsk+g/+gfy6oltin5jHa3YetXRKmOabfM8xtpmv8c8j2XkeStGPrYj2uXNWmYU9WgHaeiK9tgPQ+VJb+mWkUwtw7l5rsP6CLr0G22yrL6hHsraTsuSpp9wTMhsog+sKnHpNQ/nyluPfmtsPmXtG9tImrYoTxpB+zm2TcpFxrEfKWuI8tahLVmPZeGjbm2Ndslae5RTH3V77HxMmnLakdnYL5RdaIjtR491YQd50X77Uv+hr7VJNtqiHu2xXZzr77IizfKyQU+s33mWsoZW/do7pCfOVbZDGfUaZzamEyNjm2N6PJYNdtpm+VAO/YTIpk9I/2STkmc+1Q7tRV9s27T6J1UU9USOzoO01/rVY3tMh03UQznn8Wib5dXjGJMt6fa7eZbVtmiP44gykQ3n6lGOtByi/+Y8zmObYjssO0netmY2UU88VqdsPY9soj3kcx5D1Kce7aCcaVFmscfYIPtpuqxfm3LbJvUVupWznihvmv3u+SQ/1H+Moz5tzX6o3hjnfoh5Q8dDem2j9Ud56pHRkLzlW/LmyUgdxuYTa0eLTSy3VMe0Tf5xjo36tWlS22L5Oi4CmcBavcHLQHfRGhpEpjvwM0DO3SBycvXCMpZVjxcSflIUy6jHhcsyylI2HnvBl20b0oO8trX0aD/1RjbagbxltJE07TDPNpKHrIw5t4zy2KEt5rUWENtkmSiPXoKLQUt+VYlLyzjJ2zZssI7WxXlm7LnyxlGPdcZYOdtsXpQ3jTiXi3kex74ijT6OaeiIN2CUsc841ib5ySbroax9a1nszjbKMfaRafYfugzWrx7rsKy2k68ebLN81uM57YgcSNf/LTOUFvM5zhd1UU+uQ1k5ek5sW6K8+bbNcxl7TmybZYOvx/rNp6zjwLKkxZDlrN/yeRygWw7aRv+bhm7KZD2xztzX2kiZzEY9+pF2UVa52F592nKZjWNMWfS0/Jd0gvXTV5nVqhKX/pfHpSmrjobYqFMelI7HnNsejrHTeTQyVo9tpqzBdsY2tlhaPsctnZaJbEyLjEwbim2r3LSV8uq2fsug3yBX89QX/UE9UQ75OPZief2PMqbHNkUbKTMUrDfWM1RW+1v52qDftspMkrc88rEd6iVff8iMlCWO5aMfkpfniFyec/uRY4OMzLN+2mOeuiP33M+WtR3qox75R/vjmKJMZtNiHXUqb5prlbH2odsysY+0xTzKEWw/x7bJMnKwnPMA54Yobx36c9ZjWcrlcWRZ9NoW7ZFx7A/L2/4opx3mqc90bV9IjO22Hz1+CABj8rQHnbZfe+lb6zZPH1GPtijDuT4hK9Jsi2zQE+tXP2UNsf7sR0N6tA8d2hTbqG7izCbmIaPNMT0eywY7I2PLoJ+gHabneFo9uXw+1w4ZoS+2bVr9WV88j3oiR30b3tavnO0xHT5Rj+WIW2PdfPP0J9Lt92gL6Z5Tp/U6jsiPbDhXj31EWg7Rf3Me57FNLcZxHsry2pjZRD22P8rK1rTIJtpDfm5b1K0e7aC8aRwvVcCGWfVaTpty2+zjIduUMz/Km2a/e24cWatH/zGO+pz/opy6cpz7Iee3zmNfxXxtk1XMox4ZDclbviVvnozU1WqjdrTYqGcpY9rWmmNjHdo0qW2xfB0XgUzg0jumnLOWnPtou+8mzM023Uki53OuDichv6Iey6rHC7aWPvV4ceiPdCiLvnjswM+6hvQgr20tPU6A6lOPdiBvG7WRNO1Q3jaShy71cJ7lsUNblG9diKijJY9eghcgllmVOv+/Zbygs23YYB1+shYlZWKaX59Q3jjqsWyM1WObzYvyphHncjEvHms7afaxaeiwj5RxQedcm2RjmayHdPtGxshmG313pz5CvrYoZ4xO61eP/uOFhE+xRD3YZnl0ENSz6uzSzQvrJt02WYY464l5HucyWU+sQxl9zHNiF/UsT15LB+kxaIcXKfaHsuYjY3/KMerxOMp5bN/4ybdl0W2aZXL/UybrcVyjxz7STu0nL7NRj34U26GcepDXx63PMupxjGk7Mtl+0gzK2Veexzqty1hZ41iWY3WoM356Lxtl8xhVtqXHGz5lie3/2EYZRZZRJh5PKqMttkM50znXRyNvy9lWbbSvyFeHcuqJ4yn7ofNx5K2eKIf++H7MWF42lInp6tEe8icFy2c2k2Raedogq1aZWWxSXrvUiz7HVmYU64rl0aUeytgPQ+VJb/mROrTf+vEH81yHo2/YFm3Ka02sS/6WxZbYx5yrzzqdI8gzaCPnyjveXauM9UPKKqePk6YtypNGsP0ca4vymbH5lDVEeetwPs167HPKqdsykV9mo5xlqdvjSX2U2Wifti80tv3ocf6UcZwH3cAwj/ZYt2z0EfXYftuFbXKQVUzTFvVYvzKxbbH+7Ef6cdYT5yrboY1Rt8faYztNJ7aNMS0eK4OdLT2t/o/yHrfabh6xvj/UDu2Qkfq0KfZN1DvrsXpi/bYN/tavPus33XP1UM41yj4izfIcE7RbedLUYZ5s9CN0qMdxhNwQG+cB7aGsIfqvaTnO9sT8SfLaaNtaeiIb9Vp+6Fw95Ns2+01m5KlHO2Iax0sZHEctxrGebJPnlrGPPc9xbAt5WZ60yIdzbYqs1aP/GEd9tikyRV8r2A+tvJxmWe3K+doWbYllZDQkb9khefJlZNuMlSXWjhabWG4pj7XLurNu0ye1LcvUeRGIBGqDd/QuJYKDKcLh2AuuoRdhU8aB6qTa0qUeb1pa+tTjAuYj+lGfeqjXDYCsa0gPMupq6XESVZ96tAN526iNpGmH8raRPHSph/Msjz3aony8WUCGoI6W/KoSlz6x46dwpsfYT8W037Zhh3W06peJujxX3jjqsWyMlbPN5kV504hzuZgXj7WdNPvYNM5daJWxPs61STaWyXpIl7+LJLLZRt8zKOPIJMuj0/rVE/2HfBe4qIdjy1OGoJ5VZ5c+5SAH0m2TZYiznpjncS6T9cQ6lPHixHNi31GV5cnLOmRMnkE75KivKhv1Og5y36uLOMp5rO54c0lZdJumbbn/KZP1IGuwj2yH9pOf2ahHP4rtUE49yDsPcUywjHr0eW2nTLafNINyMs3nliP2xiymcRzti2zUGTeEZKMO7eWcvFh/PCY/suGcYP/HNspoyN5Vkqv+t3SabvF50wAAMIlJREFUn+vP6Zw7nxtbhti2OufYV+Sp235SPo6nlh8iK9eoJ7c13qjG/pHNkB7tIX9S0P5oy1D5STq1TVYtHfZxK8805Vt2ubGRGSlLrB0co0s9nNsPHBtiedJafqQO5xrrpz3m6RPG6LIt1pHXmliX/C2LfOzjqM8645gjn6CNHCuvT7pWGWsfZe3b2EfaojzlCNG3tcV6M2PzV0mu+h/lrcP5NOuxz+Gjbm2N/GyL9TkPW5aa1W1MmnLakdnYL5Rdk6A96HH+lHG03740D7u0STbaoh7tiW20PbKijGnRFtKt3w1k0gyx/uxH+rH2qCfOVbbDPPXGONsT83IbYx7HssHOlh79JrLJOjiXTSuPNMe6/pTLaYeM1KdN0+rP+vK5eiJH2wZ/61fO+k33XD2WI7aPOLY8xwTtVp40dZhHGkHboj2OI/KH2MiWMjlE/815ng/ZQ/4kedtq21p6Ihvrs/zQuXrIj3Mc55GZerSDfNM4XsowbRxZl/Vrk+fm28ee51g507M86ZGP5Ygja/XoP8ZRn22KclFfPM79EPPysW2MfRXLaFu0JeZPk7fskDz5MtKGVhu1o8XGOpY61i7n/Kxfmya1LcvUeRGIBNb6DV5/QfNHP/pR5DI+Nn3rrbcep+UDdVxwwQV9Fr96GG9KLr744u6iiy7qn2RhMHNhc/Ob3zyr6dTjp0hOpNqgHhdxLyKzbUN6qBDb+HVh7Ml6uEiPtqlHO5C3jdpImnYoTxvVjW3qacnTNtunfOtCQh25fmVh44SYb56ol0AZJ1MneduGHuuINyurJLuOdihD2/wFTeWNox5lY4ye3I/kR/lY3vbFtNaxtpOn/5mGDuo0YL+/1E2aNsnGOOuhbFwk9ZVoI7pvd7vbUXT8FdnIxP5TD+W23HLLMRN0Uj+xvH2CN+rBttNPPx3xcYh9RKJ+KQfSbBPHhmi/aTnOZbKeWIey1u858S1vecv+NMuTmHXY/l5g9E+/gY03foxl5hplo52kU9YnmdQT4yjnsX104YUXxqLd9773vc5fmcc2dDOPxTo5znrieNTXnA/jWM9s1ONcE/XQ7qF5TG6ZjWPMOvHVbH9ssPXbV57H9mpT9OeoI5aNbNQZnxiM/tsao7F+j9Vj22Ld9n9sox+WkMdfDrIj3aeEchnOc/2WMZ1zmbTq0g+ca+wr5NQR5yzSHU+kZz90Po681ZPnc39wy/EE66if45YeN6P7whP+Wa99k4tGxvE4ltM2x1jMi8e5beapVx8nXbti2+zjIT3aISP6TT3ozHMEaVE/59bBsUEdzjXWz9g0jzTa73hFVr9x/shrjfqJ4T80R5DfYtMaR9qIjBsr2uZaZawfUhZ2BPzGsSYb29sXGP3Ttzm3/dabGbs5qSxxlLcO59Osx/5ozedxrYA1QXuch+M4ULdjnfK5jzIb7aPsmgTtQY/zp4xj/22yySa9etg7jqxbNo5R9Tif63P6CH3pnIPSzEY91u/GhO1Tj/XbX/qRfpz1OFehx3bEPlK/cWRjmrFt9DzH2kbbWnr0m9j/WQfnsmnlkeZcMpSvHTJSnzZNq39Ir+nq8ZzYtsG/dV3pPKgfIRP1OJ71Q/Jth3MxcZSPOpwr9D8ZMUad6xhHpg+x0Y44Hqkn+y9prWCbWozzPKT8JDZRT2SjLH0rH8eIecTaw7FjS0ZRN3q0Q0b6DbJLGRxHmbHsrSvalPudMrbH8jnWf0hvsSE98uFcG2QtE+p3HjaOexC2SXl0xWAfkabvxfyhY+dB179czjYO9ZWMhuTVNyRPvoxgQchtbDGKbHqhZfinXTLIVZg+qW1Zps6LQCSw1m/wbr/99j2PI444InLpj5nAjzzyyP74tre97Wr5Jmy88cbdVltt1W+68Ek8F4Jf//rXze6OOuqofhFjQKOTsm5ajQuNDtTjBOQFpbapx0nTi79s25AeLgzRuc8++/T2ZD0snNE22WgH8mws0UZtxH7tQN42qhvb1NOS/8IXvjBmrLwLeItNrj+yga0LQpT1GH6U4SZN+20beuDmL2SyIMZAOyy76aabdve///37bNOM0bPRRhsN2oGe3I8oivLWG/3PtKE4MkYX/mcajKM/Xv/61x9/NRd92iQbGWc9mS2+AkvHCLpgc7Ob3az3IxnL1jHixQTlCeiQif7j4ke+m5lRD7Z9/vOfJ3scYh+RqA/KQf+PLGZhnMu09NznPvcZ2+GB9Ud/vvvd796zka1lifG/yDjKkR8ZkQc3fmE19jWMYoC5F7wx3WPZRLb2v59mW/aQQw4Zb9qR1ur/lh7lifU1Lhjp49jGzCbPY15k0h7aPTSPUU+LDcwJXOAS8NXsv33G3/5Zv30VWcVyHKs7pme/iWzUeetb33osEv0X23yqiALkxfqzbW7cjpX97SD3ERteBC6Y4zj4W/F57chzoGWIc/3maaPnxHCgrtjX+gF5sa8or243lEgjOJ6cx1elrvJD5+Po/+qhjhiOP/74/tTxpF7XMTKjHtuk30RdWTd51msfx/IcR1+JTGI5baP/WtcKlm3VH+vAj5S3HbFtbGxQx5Ae7ZBR9EPmwTxHoCeuB9jS8iMZOddY/y9/+cuen+swa0VkpN8wF2B3XmvUw9iB86Q5osXGvnHO5CkbbLQNbAARsCmuVdSLPfohbPRfbMX/Ixvk41wvX3RnNplxa4NX+ViH86mM1eMTQa353HZjB6wJ2uN6LiP4oTuv57GPZBTZRP/rKxj9s988nxRHP95ggw36ojKO8+Dmm28+VkP9ca6XDW2hb9SjrznW9ZE852Q26vGpYTeINYBrrli/cw31RT/Oepyr1EPZ2EemG0c2phnHtcY04+g3tK2lR7+RjbLEsf9kE/Pjsb4q65gX7ZCR+rRJ/8ty8XzSsXqizc6D6G5dVzoPwj/Pp9RlW/TD2I5oS5Qn3bHlXGFZ9BAYoxw7jmSX2ajHNimvvuy/pudYNpkx48D5I8tMYhP1RDbqoG8t41gzj1h7OHa+lHX0Q/Roh4z0G2SXMmy77ba9Ohl7nXbqqaeuNg60Kfc7CpwrWrZl/2mxQY5+Z60xyEbW1q//GGd7nP/0H/UZ20ecy9e8SbGbwXk+RCa2caivZNSSt63oGpInz7HhGMttzIwyG3QsR9C3h9ZDr6MmtW057CqdVxwC83exrjjtmrkld77znbvNNtusO+2007rDDjtsntxHPvKR/gkVnsrYbrvt5uXlkx133LFPcvL44Ac/2G8k/PznP+/222+/Ps9PHB/5yEdm8fG5erjgYBHnBgPb9t9//7EeLrS5yGeBG7Ktpcf6Dz300L6+qMfFM9omG+1Q3jZqNHYobxl0a5t6zFOeNvJpOZ8MZ3l05wtZ2xTlM5u42KkTXbEfHvKQh5DUf2pv29QTFwEv5NDJzcKvfvWrXo7NPNukvDF6dt999/EGci/wt3/ctKJHf4ARCzecorx+qP95cxnbE/VyrD0yxv+22Wab3rdlrAz1unjYNm2SjYzVw9MwyMWAr2ij6W502lf0MUxok2nxYgG5L33pS70ejrXfGP/PbNSD38lGtvYRugwtNl4cZvuVibFl7Cttc4xT9owzzogi8461CRuZR7Q/yiOAj8aLFpVwsRb91/rtK/vIeezoo49WtKOPrL/lP7LJfWT/owgf4UMdn9ZSeex/2UzSg9yPf/zjsf/bDtKpo8VGVlxEa5OMJs1jLTb4ETepXpjjq/atfpQZWf8kxvioTxhEeXUPsUGnG0aysa+w7Stf+QrJYzZDfYWeyNJ5C9nYR64Z2EtAznHAOT5m/ZxPC5GNerDRm4UoT108iURwrOY5J9qj7qgjjwPzaGNmY15Lz9lnnz24nua+Qo+6W31sHw2tVbFN6JqFcRzrtG2W4NjIdTgfk247GKMnnHDCWG2LkXrso0nrOWUdU/r8WPmEA+uNtnMjx7WO67DfykFNXqtkox7WGtco+6o1R2hSi41rjXOMthljnze4jjnj6IemWRd+sPfee8+71lGnZWJsm2IZ/f9rX/taLDrvWP5ee7X0sAnKet6azx3/rlUqVw+M7RtvUON6nvtIDpGNa4S6iYfGUSzj8ZAfw9j6KPuzn/1MkXnzYGaDXNxAYK3Rn/GROB5R6Fqlctmgh01cAte0bohyjl/a/87DrhXZj6Me5qrFrOfUbTj//PM9XC3OftNizNe4ucZwbGk/yuy/zGa1ikYJrrUypoxrg3bIyL6ijDbpf6QZrN/NNtNbsXqcKyhD/fZj67rSedC+QgY9bNznQP/Z15ER5aK8ctbrOTFtxH8co/o1NrfYIBP1OH+Rnv03r1WUMcgmM3asUy7KR92xbUN6YJPvPZxz4nwc7ZFxnGeYh/TDfD0Bo+g36lqqmG8cEmTsJiRrTW5by2+0I16nmeaYcByY3mJjHv7QCu9///sH76tiX7VkYxp9bB/F9FmOmQ8IcawpZxsn9ZWMWvLalNcq9cc4jg3TXU/tI8fYQtioa01ix0i+HkPXLGzWpM6SWbsI/MMeo7B2NXl+a1msuOH98pe/3F+ksAHBRRUD7NOf/nT/FNWrX/3q/qnM+ZLzz7gx5EmAM888s5dhYf7Upz7VHXTQQR1PpXDhwVNnd7vb3bonPelJ8xbJqEk9XpzyhAThG9/4xrybcT5NReeQbUN60BVvGvxUlrRsW2SD7S7scaFFHwF5bny0lwub17zmNT23qIdFEXku8uKkneXRySLOJLzDDjv0emyTjF3gZUP9TtLIUy/1ofvggw/uL3Ro46677trzbDGWB/Lq5wX0n/3sZ/v+I33nnXfuPxXUb2Rj2/GhVqC/DjzwwPHFrIzYiMX/1PPVr361t5cLe9uAPi74vBHhmAWOVy086EEP6stpD23mBuMzn/lMz1gZbeJcG20bi3hkI2P8mHHARa8XH+r54Q9/2NvtOTFssMm+kjFtYqOZvqc/YqAd+gS2xbGy00479eMysjn55JPHHLzI4qKLMUs5wyQ2+AN/XxptLg8FfPSb3/zmeB7An+0rGauH/stBeRZwwi677NK/okE2P/jBD/qnDugrbtrf/va3z7NfffQB9eR5RD9GD+3gw4Of/OQn/YY6slxcwcP6o/+02HiThIz9hh76i3ER+1/dzpGRzZAedMETX4pzBek8ifeIRzxi7DeywT/POeecJpfoR45b0qIfZzZxLLApf9xxxzXHWJ5zsh6ZwgGW9CGBMQ0n+uFzn/tc78tDbNjcdRwiCxvbxA0W6wlBNnEejYzRY/sp77zF8SmnnNJvECBLGcbW61//+r5eGQ/5n3OMMW1jzm6x4SkI9GDzd7/73XlzO3Xiw9roWM1zTrTnvPPO677//e+Pb6JoSxwH+A9jkE2Hl73sZf164vwHG7iyjh9zzDF9jHwMrhmmySb3FXrwP3S1+tjNtqG1KrapNcaVc2xhx1577TVeq1rXClxXMBcQlGfM4jv8MVc4DzofUzb6D31EYI6CK+2L4Vvf+la/0erGi2sV306Jepx/eTKcD/Fa86n+45xjPc6DrDX2J3n0jb6iftLzWiUb9cQ5i/IExxPH6uSYMMRGdpSJ8wXnBGyCW16rsId0/TDaju+cdNJJvTw2MSbscxJnYdPy/17h6J9rjfOx14WyiYzRg+2xfudz5zXXKvWrJzJmXMA0MuaY6468VmU2jFHmOOYq5so8juhrPzTKbKL/6cd5jsFuv2pv24bYMEaZJ2mLc5ztZlwccMAB42s20p2PLSMb9MTrm+hv8Lb/TYdVvOYe0pPbYR+5njsfY09kM8RYu2HPfB7XKv0m6pExPsaGy6R5MLOxLuKzzjqr3wDLjMnja9LMW17zwog+1R7KaJObaNF+/YfNHvrAOiIbdBCyHtKof9J1pWuVc456mB/xnRjyWHduc62i/hjod75lp7+bp59wno8zG8q09OR5nHKuGbOwkTFjlXkjBmwaWqtajB2jrKkE53Pn2zgf9wVG/9AjY++v0INdBs5b91WsVTHof65H+kica8zLc07Uw3jRB01nHPFhZG5by2+UYT5yHTdtaB1tsVGGe9Osh7xJ91XRjykrG33DtcE+to8oO4mNemSLz2e/Hlqr0J1Di5HyzoN5rco6OG+Njbie4kO0Mc7HLT0xzeuw6D+T2ERZjh0j3DvGOTbec8b5L8tPOrcfFmLPJH2Vd/kksNZv8NJtfIJ0m9vcpjvxxBO773znOx1fi+LCbrPRp9U8jTnLI/JMOnw9j8WdjQMu3lwcdY1/+7d/617wghf0Fy6m5Tjq4aI8XjzHstNsm1UPOik7ZFtkw4XTpOANDZPpHqPPDSI39fDkEBOpZaM+0/jElsWcSTduKMQ2ybglzw0W8uqjPi5QbCMX/PbVJMbqZpFDnsWKDeQHP/jB/YazbcJvprFBF3q8MImMWnrijSUc8AMuPlwsWxci6pExNkc9tifGtq3FRsbY7IIfZW0z/ZLZxL6SMRdA9knU43HM054ttthiPDatr9Umnk6JfYTOuLi12Mgy2s+vWauf/DwPTNPDDZB9rLz99+hHP7r3m8iGm4Nvf/vb/QKf7Y9cGAcG2UQ/Rg+cY6DPsN/6o/+02Dj/eSGddeXzITaT9KjDvnasw5UPKjIbborpj3xTpJ4YIzsrG+TgbciMvPnJ9rQYM+dzQcsYdZxw8RjnaP1mFjbYhDxfp2Pukg3pC9WjLyOrPbzv0flvkv85xxh70T7ERj+GA/3rjW1ewxyrQ32FHt49SJ/rJ9jPcRwHpPF1OuZjQmbjOq6P9YUG/smGNaulZ6iP7SP8Z2itmsT4Jje5ybynY7hBxC9lg94cvLEgXXkZu9bl+Vgdto0PNpj7naMyI+YA57G4VmU9xx57bL8m0jet+ZTy+k+cc0iPY4tNrNjX5OeQ1yrZRD3e9GVZzyk7jc0s63m0NfeVjGVj3TnWDtJnYRP9n7HRWmuiH6M3spExepynKEPIa4Vr1arc+XpkHK+xLEccx+gQGxi7KRHnOsfRtJtmGevH+n+0w2PbNokNN9hwQU9ca+I4cIxQN2uVITJmk3eWgExmM02P7YARdrqeOx9br2yGGLuJYz/ltSrrkTHje9o8mNmoi9hNB9eReK1FfmutivcQlEH/4Ycf3s810X7nYfqIjS3ryGzQQYh6VqXMr9+0vFY555iPzWz2udaZ3orjWhXz6Xf6svXhWCzncfZj01t6Wv7rGJuFTWSs/1vftLUqM85j1PncedD7KvUby9hrrJaeSWuVevQ/1yN9JM415uX5WB3E6snjKJaxbXmMxzKu43EcMAe11tEhNuhrbYLGetBpGLLHNnk9YXn72D4ifRY2mW1so9ccQ35s3cQtRso7D+a1Ksp73BobcNG/KTfERh051rboP5PYZHnOGSPsPcU5Nt9XteSmpdmfC7Vnmt7Kv3wRuNLIyS8d/Zcv25fFWi5geCKKCwY+AWRiWGhg0mBzjJiFh69tMXHG94PNolM9XLCwqPGuMxYO7FqIbS093AiyWcnXaZlkZrFNNixsymMTF57owkZ0TbMNPSwcLNxrIg872yTjVv2WwcWH+FumxZiNYjZ0Jsnbj5ENfUSf21fqsY5JjLKehfZ1tEfG9jH26NvaNKltspGxfoyerHuS/6gntt/61c1iPW2stNioZ1I7ZGIc/U822X7rmjQPTNIzizz2yKZlv3kykn+2taUn9vW08SgX4mi3bHP9s+iepKc1VqMNHtt+2bChzbyqHzF3OA8N9SO6sh7tZx7jeNYxNqQn+khs9xD3WEbGtlHbok555HgWPdPamNvU8q1cb+t8SM8kG1t1tfSYRjxtjtC2WK/9oJ7oP9P6v6WnlWa9rdh67eNZ293SNZQ2Sx1ZttUO9cho0lqlPvQsZD1QLsbWS0wfxzkij5EWP3WpR/tnnSOUN5ZN61pnlrUq6pGN1zqT5nHlYmybZJPltXXanDGkZyFzDnapR8b4iH00ZGNsj8ez2m35oTjq0Q7H2qxts03K5bVmlnGAfVmP9Uc/mrRW2cYhPdP6WHnjyMZ5cJY8yxi39LTSLL+QWD3T1qqsU7mFMhnSE+vPfjRpzlEf9iz1WHceZO6YtlZpB7H+E8do7v9Yfui4xVjdjpVZ2Szk3mMWe+yjxbZxqK5Z0yMjbVoIG+uJ/jPLHKFcjO0bYtbRhaxVLT1r0o6oJx9HVkvlj7mOaeeZUV5Pp8kvV/5i2SyXXaX38kugNngvv31XlheBIlAEikARKAJFoAgUgSJQBIpAESgCRaAIFIEisJYTWPjjqWs5sGp+ESgCRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFYKUQqA3eldITZUcRKAJFoAgUgSJQBIpAESgCRaAIFIEiUASKQBEoAkVggQRqg3eBwKp4ESgCRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFYKUQqA3eldITZUcRKAJFoAgUgSJQBIpAESgCRaAIFIEiUASKQBEoAkVggQRqg3eBwKp4ESgCRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFYKUQqA3eldITZUcRKAJFoAgUgSJQBIpAESgCRaAIFIEiUASKQBEoAkVggQRqg3eBwKp4ESgCRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFYKUQqA3eldITZUcRKAJFoAgUgSJQBIpAESgCRaAIFIEiUASKQBEoAkVggQRqg3eBwKp4ESgCRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFYKUQqA3eldITZUcRKAJFoAgUgSJQBIpAESgCRaAIFIEiUASKQBEoAkVggQRqg3eBwKp4ESgCRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFYKUQqA3eldITZUcRKAJFoAgUgSIwE4Hf/OY33cknnzxT2cui0LHHHrtatZdcckl3pStdqf976Utfulp+JSyewHL5xUrou5vd7Ga97zz84Q9fPKjSUASKQBEoAkWgCBSBInCFI1AbvFe4Lq0GFYEiUASKQBG44hI44IADuq222qo79NBDV1wjf/azn3WPecxjuu23337F2XZFN2gl+8UVnX21rwgUgSJQBIpAESgCReCyJ3Dly96EsqAIFIEiUASKQBEoAtMJsIH6qEc9anrBy6jEHnvs0e27777d1a52tdUs4Ondq1/96n36la9cl1+rAVpEwnL7RfXdIjqnRItAESgCRaAIFIEiUAT+LgTqDuPvgrkqKQJFoAgUgSJQBNZmAuuuu273+9//fm1GcLlte/Xd5bbryvAiUASKQBEoAkWgCKw1BOoVDWtNV1dDi0ARKAJFoAgUgSJQBIpAESgCRaAIFIEiUASKQBG4ohGoJ3ivaD1a7SkCRaAIFIEicAUk8M1vfrO78MILxy07++yzu2984xv9+R3ucIfuH/7hH8Z5Hlx88cXdSSed1P9xvO2223a3ve1tu0033dQig/Fpp53W6z/jjDP6Vyvc8pa37PjbfPPN+x+7ioLnnntuhz0///nP++S5ubmxbRtttFFf31//+tfuW9/6Vp+/8cYbd5tssslYxU9+8pPu/PPP79Zbb73+/cJkoOu4447rvv3tb3c3velNu9vf/vbdrW51q26ddaZ/Nv/1r3+9O/7443ud22yzTXeXu9yltwG74EjYYostug022KA/Xsp/f/rTn7qjjjqqgx+vTuDHwWR37Wtfe2pVC+2zNfGLqUakAn+vvqPtn/3sZ7vTTz+9u971rtfd/e53725961vP1OeavFB+3/ve9zp+nI6w2WabdTe4wQ1UNS/+8Y9/PPbvm9zkJt2NbnSjefl1UgSKQBEoAkWgCBSBInAZExhd7FcoAkWgCBSBIlAEisCKJjB6b+3c6JKp+ferX/1qnu2jDbm5t7zlLXOjd+E2y//Xf/3X3C9/+ct5Mp784he/mHvYwx42N3rvalN29ANqcz/96U8t3sevec1rmmWx97nPfW5f5te//vW4zEte8pJ58s95znP6vHvd615zow3SuZ122mlcNrZ5tOE3N9pInicbT0499dS50WZuU/apT33q3Ggjb5z30Y9+NIouyfGHP/zhudHm9biOaPvoNQdzH/jABwbrWdM+W4hfDFY+JePv0Xf40GgDfDV2ow8len8bfbDQ5+GbrbCm/D70oQ+N6xx9iDD35z//eTX1+Nzow4C+3IYbbjjRB1cTroQiUASKQBEoAkWgCBSBvwuBeoJ3dPdRoQgUgSJQBIpAEVjZBHhK9w9/+EP/NC6W3vjGNx4/RRif3uU9tw95yEO6I488sm/QDW94w2677bbrrnOd6/RPtfLE4mhzszv66KO7z3/+890tbnGLccPR/9CHPrR/cpZE5Lbaaqvukksu6b74xS92PB1JzJPAJ5xwQv/EI+V4mvGOd7xjd9ZZZ3UXXHBB/4Qv9hJ42nHW8Nvf/rZ7+MMf3h122GEdT7vyxO5ow62vi6dIsflOd7pT96Mf/Wj8g23qPuWUU7p//ud/7kYb1H3Slltu2dvE08U80bvnnnt23/nOdyy+5PEnP/nJ7vGPf3yHneuvv35vC08I85Q1T1HzlOguu+zSffe73+3e9ra3zat/MX02q1/Mq3AZThbTd09+8pO70eZ3bxXv++WJa9jxBPfJJ5/c++Hvfve7QasXw++xj31sd8ghh3Qf//jHez9785vf3I0+gBjXRX9S5qKLLurTPvjBD/Zjb1ygDopAESgCRaAIFIEiUARWBoG/yzZyVVIEikARKAJFoAgUgUUSGL3GYPy04etf//qmtle+8pXjMk972tP6p1Zjwf3222/uWte6Vl/mfve7X8ya+9jHPjaWfc973jMvjycb3/jGN47zd9ttt3n5nPCU7Ojqrn9yOGfO8hQosvw961nPmqO8gaeKebrS/L333tusPubpzdHGb59/latcZe6d73znvPwzzzxznjx6PvKRj8wrs9iT0UZ4X//olQyrMT/nnHPmRq+Z6PNhP9own1fdYvoMRbP4xbwKF3iynH136KGHjvuVPoxPaNOvL3vZy8b59FvrCd7F8htt3s6NPjDp6+Gpd54EN0Sf92l08youAkWgCBSBIlAEikARWDkEpr/IbXQ1WaEIFIEiUASKQBEoAiudwOjVCd1oQ6o380EPelA32qTtRhuK88z+z//8z/5pVhK/8IUvdJ/+9KfH+V/96lf7Y56e5KnKGEavAuie//zn90/88p7SH/7whzF7yY7vfe97d+94xzv6J3hVyvtYDzroIE+7gw8+eHzMAU/8+m7dF7zgBd0zn/nMefm8W/Xwww+fp3NegUWe8NTw97///V7LaHN6Nea8c3j0CoJutHnYv3uY9wMbFttn6lkJ8Zr03ete97re9Ktf/erdZz7zmXlPx45eE9JzG72yY7B5S8GPJ65Hr2ronzz/4x//2D3pSU/qRrcq3Yknnti9/OUv7+vmHdCOrUFjKqMIFIEiUASKQBEoAkXgMiNQG7yXGfqquAgUgSJQBIpAEVhKAmxS+VX2V73qVYOqH/3oR3e8uoGAjMEfjuLr6AceeKDJ45hNXl4xwI+HjZ4EHqcv5cELX/jCpjo2afkRNgKvgYiBVzAQ2Mwekh+9O7Ubves3ii3ZMXa5kQ43v84fK4A5fcNG8Ohdw+OsxfbZWNEKOBhiP9R3vP6ATVTCzjvvPPbJ3JSXvvSlq/2wn2WWit/o3dLd8573vF7tMccc042eAu/YWOYVIbwuZP/99++uetWrWm3FRaAIFIEiUASKQBEoAiuMQG3wrrAOKXOKQBEoAkWgCBSBNSNw2mmn9YI8DTl6JUA3+iG15t/oR9n699tSWBmOeerXwJO+d7vb3brRqyD6d5PyRCMhvu/XsksZj15xMKiOjUICm24x8D5gAu8L5mnMoXDPe95zKGtR6Wz88f5fApuDvP/3CU94QnfAAQf0/EnnadR11ln9slP+a9pn6F4pYaF9d8YZZ3S8u5dw17vedbAZW2+9dcdT0K2wlPx4mvjWt751Xw0fBvC+asJ73/vebosttuiP618RKAJFoAgUgSJQBIrAyiRQP7K2MvulrCoCRaAIFIEiUAQWSMDXJvBjabxGYZbAJhubt2xA8jV0fuxq9O7e7v/9v//XHXvssf0fT1ButNFG3QMf+MDuUY96VHff+96342ne5QibbrrpoNprXOMafR5PfsbgBq8bwDEvHm+++ebxdEmP99prr/4H6vhRtQsvvLDjx7j4Y0P8zne+czd6d2y34447rvajc4vtsyVtxCKVLbTv+PE5wyabbOJhM+bH+kbvMl4tbyn58QqN0buZ+x/y41UNhMc97nEdT19XKAJFoAgUgSJQBIpAEVjZBFZ/lGJl21vWFYEiUASKQBEoAkWgSWD0A1XN9EmJv//97zveIWvg/aNf+tKXukc+8pEdT5UaRj/k1e2zzz79U773uMc9unPPPdesJY1jnbMoZiP6vPPO64u6ATwkt+666w5lLTqdDfWjjjqqe+1rX9u/p1iFf/nLX/qnenfdddf+KdA999zTrD5eij6bp/AyPFlo3/HhgoH3Pk8KfMDQCkvN7/rXv34X/WT0g3itaiutCBSBIlAEikARKAJFYIURqA3eFdYhZU4RKAJFoAgUgSKwZgR8CpKvtP/617+e+Y/308bAqxl4vQDvkv3c5z7Xv5sUnYavfe1r/VOO+VUJ5v89Y54k5ulOAj+4NSlMy58kO0veNa95zY6nnXnP7umnn969613v6p/q9f28f/rTn/qno/nxO8NS9Zn6Lk8xP55nmPaBAU9Ft8JS8+PVGrGuT3ziE93ee+/dqrrSikARKAJFoAgUgSJQBFYQgdrgXUGdUaYUgSJQBIpAESgCa07g5je/eS/Me0mvcpWr9D8OxQ9ETfvj9QytwBOxD3jAA7q3vvWt/ftI+YG1W97yln1RNuR4hcNKCLb7zDPPnGjOj370o4n5S5n5j//4j90znvGM7uCDD+6fkOZdxoaPf/zjHnbavlR9NlZ8OThwYx5TzzrrrIkWD+UvJb93v/vd3WGHHdbbwUb9tttu2x/zPt6/p+9MBFGZRaAIFIEiUASKQBEoAk0CtcHbxFKJRaAIFIEiUASKwEojEDdi/dGzaCM/MkbgtQWHHnpozJp3zDts73Of+3T3ute9OjavCLxKYOedd+7+6Z/+qXvsYx87r7wnbO6+8Y1v9LT7+te/Pj7mQPtats0ruMQntIXwk5/8ZLxB16rife97Xyt50Wmf+cxnuh122KHjR8ZOPvnk1fTx6oIXv/jF3R3veMc+j/f0ymgxfWZFcudcveat5JhXffDUM+HAAw8cNJUfOxt6+nop+FExT12/8IUv7G247W1v2+2xxx79k7u8Q/k3v/lNt9NOO/VjZNDIyigCRaAIFIEiUASKQBG4TAnUBu9lir8qLwJFoAgUgSJQBGYlEN9xGr9GrvwTn/jEbr311utPn//853e//e1vzZoX85Vz3rP7la98pbvqVa/a57GR9eMf/7j75je/2e27777dt7/97XkynsQNxLve9a4m97H28SoCNsX+XuHZz35259f9X/KSl8z7ir028ONZvCN3OQLMefKTJ4h33333wSpkd5e73GW8Gb6YPrMiuXPe8gvLrbSYzV1+tI9wxBFHdEceeWTTxN12221w43op+PGqEX5IjfdR8+Q7P45HfIc73KF70Yte1Nt03HHHda973eua9lViESgCRaAIFIEiUASKwGVPoDZ4L/s+KAuKQBEoAkWgCBSBGQjwLlc2ngi8G/SQQw7pNy3ZUCXc8IY37F75ylf2x3ylnadxeVrUwNfM3/SmN3VPf/rT+6TrXve63bOe9Syzu8c//vH9MRuRT33qU7tTTz11nMfB0Ucf3b92gGN+jCpv8K6//vpk9QE7vvrVr/avdjBtuWJeQfGKV7yiV3/SSSd1t7/97Tue1j3++OP7jdenPOUp/ROYsf741GtMX5NjnkTllQyET33qU/27d+P7iXkf8nOf+9zeHso84hGPIOrDYvsMJdP8YlVNK/P/O9/5zu52t7tdb9zDHvawbv/99x8/KXvxxRd3j3vc47qDDjpo0Pil4McG8gknnNDXwasZeILXwIa9759+1ateNW88WabiIlAEikARKAJFoAgUgRVAYHQTU6EIFIEiUASKQBEoApcLAttvv/3c6PJp3t9oc2ps+2hjcW6XXXaZlz/aeJ0bvat0Xtro6cm5Y445ZiznwWiTd1xutAk6d+Mb33hutAE3N9pIG6ff6EY3mhv90Joi4/jLX/7yHDLRvtGmXZ8/2uQcp4+esh3LcDB6TcQ4b/SqiHl58WS77bbry41edRCTx8ejHzWbG/3o2lhXtGP0Y1xzo03ncd5og3wstxQHp5xyytxoo3msf7TpOrfNNtvM3eIWt5gbPWE7Th9tNs+NNuTnVbnYPkPZNL+YV+ECT5a7784+++y5jTfeeMwIjqNXL8yts846fdqd7nSnudEP//XH+lNswmL4jZ5iH9dzm9vcZrW+oR58XVsYR6On02P1dVwEikARKAJFoAgUgSKwAgjUE7wrYJO9TCgCRaAIFIEiUARmI7Dffvv178711QpIxSdtRxuc3fvf//7u8MMP70YbjN1oY6r75S9/2fEjXgRexcD7RE888cTVnsAlH9n3vve93WabbdZ/Lf6cc87py/7sZz/reOL3oQ99aPetb32rG222UnxeuOc979k/vTrarBunR9vGict0wI+a8RqGpz3tad1oU7D/cTmeYuY9wzzNy5O9hutc5zoeLkl8q1vdqhttBHY77rhjz5zXY/CjdD/4wQ86nrDmKdA999yz//MpbCtebJ+hZ5pfWNdKjEcfInSf//znuwc+8IG9eZdcckn/Tly44G9f/OIXuw022GDQ9DXlx5PVj3nMYzreSU2f7LPPPuMn5GNl+DqvPCEwjjyOZeq4CBSBIlAEikARKAJF4LIlcCU2mS9bE6r2IlAEikARKAJFoAgsjACbhmeccUb/qoQNN9xwUPh3v/tdvwF87rnndptuumn/Q2C80mBaYNNr9GRl/8Nlv/rVr/rN4s0333ya2Dj/vPPO66h79ORsd7WrXW2cflke8B5eNrcJecN3Ke1i45BXZPCjb7y2YvRkaP8ahVnrWNM+Q/+sfjGrLX/vcnygwA+e8QEGm/Tx/cKz2rIYfrPWUeWKQBEoAkWgCBSBIlAEVhaB2uBdWf1R1hSBIlAEikARKAJFYEEEDjzwwP7J3S233LJ/ejc+3RwVPfnJT+4+8IEP9E9pXnTRRd26664bs+u4CBSBIlAEikARKAJFoAgUgcspgStfTu0us4tAESgCRaAIFIEiUARGBHgNxXve856eBV+190fkIhx+8G2vvfbqk+5yl7vU5m6EU8dFoAgUgSJQBIpAESgCReByTqA2eC/nHVjmF4EiUASKQBEoAms3ATZseafuxRdf3O2+++7dBRdc0N3//vfvv+LP6yWOOOKIbtddd+3fKTz6Ebhut91264GNftCtG/1A16Lg8U7j/E7dRSlcQuE//vGPfZsXo5LXa8CsQhEoAkWgCBSBIlAEikARWMkE6hUNK7l3yrYiUASKQBEoAkWgCMxA4NBDD+1/kIt3Bxt4sjeesxH7tre9rePH2AgcL/YHsx796Ed3++67r1WuqJhNb94HvJjAj4ptscUWi1FRskWgCBSBIlAEikARKAJFYNkJ1BO8y464KigCRaAIFIEiUASKwPIS2GGHHbrvfe973Rve8Ibu4IMP7njHrpu76623Xvev//qv3bOf/ezudre73diQ613vev2Px40T1uCAH5FbqWHrrbfufvOb3yzKvKH3GS9KaQkXgSJQBIpAESgCRaAIFIElJlBP8C4x0FJXBIpAESgCRaAIFIHLmsAll1zSnX/++d2GG27Yrb/++pe1OVV/ESgCRaAIFIEiUASKQBEoAstIoDZ4lxFuqS4CRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFoAgsJ4F1llN56S4CRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFoAgUgeUjUBu8y8e2NBeBIlAEikARKAJFoAgUgSJQBIpAESgCRaAIFIEiUASWlUBt8C4r3lJeBIpAESgCRaAIFIEiUASKQBEoAkWgCBSBIlAEikARWD4CtcG7fGxLcxEoAkWgCBSBIlAEikARKAJFoAgUgSJQBIpAESgCRWBZCdQG77LiLeVFoAgUgSJQBIpAESgCRaAIFIEiUASKQBEoAkWgCBSB5SNQG7zLx7Y0F4EiUASKQBEoAkWgCBSBIlAEikARKAJFoAgUgSJQBJaVQG3wLiveUl4EikARKAJFoAgUgSJQBIpAESgCRaAIFIEiUASKQBFYPgK1wbt8bEtzESgCRaAIFIEiUASKQBEoAkWgCBSBIlAEikARKAJFYFkJ1AbvsuIt5UWgCBSBIlAEikARKAJFoAgUgSJQBIpAESgCRaAIFIHlI1AbvMvHtjQXgSJQBIpAESgCRaAIFIEiUASKQBEoAkWgCBSBIlAElpVAbfAuK95SXgSKQBEoAkWgCBSBIlAEikARKAJFoAgUgSJQBIpAEVg+ArXBu3xsS3MRKAJFoAgUgSJQBIpAESgCRaAIFIEiUASKQBEoAkVgWQn8f3o9ZmRat/ycAAAAAElFTkSuQmCC" />

<!-- rnb-plot-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5nZ3Bsb3QoTUNDX3NkX3F1YW50LCBhZXMoeD1uYW1lLCB5PXZhbHVlKSkgK1xuICBnZW9tX2JveHBsb3QoYWxwaGE9MC41LCBjb2xvcj1wYWwoNClbMV0pICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgYXhpcy50aWNrcy55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCIsXG4gICAgYXhpcy50ZXh0Lnk9ZWxlbWVudF90ZXh0KHNpemU9MTQpLFxuICAgIGF4aXMudGV4dC54PWVsZW1lbnRfdGV4dChzaXplPTE0LCBhbmdsZSA9IDApLFxuICAgIGxlZ2VuZC50ZXh0PWVsZW1lbnRfdGV4dChzaXplPTEyKSxcbiAgICBheGlzLnRpdGxlPWVsZW1lbnRfdGV4dChzaXplPTE2KSxcbiAgKSArXG4gIHhsYWIoXCJUZXN0aW5nIFNldHNcIikgK1xuICB5bGFiKFwiU3RhbmRhcmQgRGV2aWF0aW9uXCIpIFxuYGBgIn0= -->

```r

ggplot(MCC_sd_quant, aes(x=name, y=value)) +
  geom_boxplot(alpha=0.5, color=pal(4)[1]) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 0),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Testing Sets") +
  ylab("Standard Deviation") 
```

<!-- rnb-source-end -->

<!-- rnb-plot-begin eyJoZWlnaHQiOjQzMi42MzI5LCJ3aWR0aCI6NzAwLCJzaXplX2JlaGF2aW9yIjowLCJjb25kaXRpb25zIjpbXX0= -->

<img src="data:image/png;base64,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J8cVV1yRue7Awje/+c348z//88xn3759B1ZJ5/8d67avueaaTDvNzc3x3//7f8/0OXMiVejq6orbbrstkiBw/5YErI/2pXH99/BNgAABAgQIECBAgAABAgQIECBAgMCxFTj6JajHtp9Za+2qq66Khx9+OPpfnvaFL3whHnvssTjzzDPTAc8kN++Pf/zj2LNnT6YPN91006AXtGVO/K6Q5LV94oknMoc7Ozsz5YGFsW77nHPOSefg/cY3vpFuJsnD+0d/9Edx2WWXxQknnBBlZWXx4osvxs9+9rNYv359pivHH398DAwOZ04oECBAgAABAgQIECBAgAABAgQIECAwoQUKPsCbrMZNVu5+4hOfSAd2k9lKArzJZ6jtuuuuiwsvvHCoU0d8LBttv+9974skt+69994bfX190dbWFv/6r/86bN9OPvnk+Lu/+7uYMmXKsHWcIECAAAECBAgQIECAAAECBAgQIEBgYgoUdIqG/imprq6Oz3zmM+lVrLW1tf2HB32fcsop8Y//+I/xh3/4h4OOj3ZnrNtOXrR2ww03xB133BEjvTRt0aJF8Rd/8Rdx++23j7gaebTjcz0BAgQIECBAgAABAgQIECBAgAABAtkTKEqt8uzL3u1z88579+5NpzJI0jLMnz8/kmDovHnzjslgxrrt1tbW2Lx5c/qT5OWdMWNGLF26NP05JgPSyLgLJPOepCBZsWJFOk3HuHdIBwiMIJC80HLbtm1x3HHHZV6COUJ1pwiMq0BDQ0Ns3749klRHxcXF49oXjRM4lMD+/ftj586d6bRdyYIAG4GJLJC8AyV5YfTAl0hP5P7qW2ELJO/cSWIHq1atKmwIo88JgeRZTf5OkPz91ZZfAgWfomGo6Zw5c2Ykn/HYxrrtysrK9H9o/MdmPGZTmwQIECBAgAABAgQIECBAgAABAgSyK+DH99n1dXcCBAgQIECAAAECBAgQIECAAAECBAhkTUCAN2u0bkyAAAECBAgQIECAAAECBAgQIECAAIHsCgjwZtfX3QkQIECAAAECBAgQIECAAAECBAgQIJA1AQHerNG6MQECBAgQIECAAAECBAgQIECAAAECBLIrIMCbXV93J0CAAAECBAgQIECAAAECBAgQIECAQNYEBHizRuvGBAgQIECAAAECBAgQIECAAAECBAgQyK6AAG92fd2dAAECBAgQIECAAAECBAgQIECAAAECWRMQ4M0arRsTIECAAAECBAgQIECAAAECBAgQIEAguwICvNn1dXcCBAgQIECAAAECBAgQIECAAAECBAhkTUCAN2u0bkyAAAECBAgQIECAAAECBAgQIECAAIHsCpRk9/buToDAeAl0dffEhm11sXHbntixc3fsaimNxfNqYum8GTFpUtF4dUu7BAgQIECAAAECBAgQIECAAAECYyggwDuGmG5FYKIIbNlVHz/99YbYvrsh6hpaoqm5JV7Z2RYzpk2JxXNr4vwzVkTttMqJ0l39IECAAAECBAgQIECAAAECBAgQOEoBAd6jhHMZgYkqkAR3v7/2xXhl896onFwaM6ZWREVJb1RWVcS+htbYW58K+LZ1xGXnrkmdE+SdqPOoXwQIECBAgAABAgQIECBAgACBwxGQg/dwlNQhkCMCSVqGZOXuq1v2xpwZVbEklZKhekp5lJcVx7SqybFiYW1UVZTHq6ng76OpejYCBAgQIECAAAECBAgQIECAAIHcFhDgze3503sCgwSSnLvb9zRGRXlq5e4wKRjmzqyO7t7e2LRjf+zc1zToejsECBAgQIAAAQIECBAgQIAAAQK5JSDAm1vzpbcERhTYVdccjc3tUZNKyzDSVlNdEY0tHbFLgHckJucIECBAgAABAgQIECBAgAABAhNeQIB3wk+RDhI4fIH2zq7o7umN0pLiES9Kzif1Orp6RqznJAECBAgQIECAAAECBAgQIECAwMQWEOCd2POjdwSOSKByclmUlRZHR2f3iNclgd3SkknpVA4jVnSSAAECBAgQIECAAAECBAgQIEBgQgsI8E7o6dE5AkcmMH/W1PTL1PY1tA57YV9fX9Slzk9PpWmYn8rHayNAgAABAgQIECBAgAABAgQIEMhdAQHe3J07PSdwkMCSuTWxOPVJgrg79h78ArXk+Oad9TGloixWLJwZtdOnHHQPBwgQIECAAAECBAgQIECAAAECBHJHoCR3uqqnBAgcSmDSpKI4//Tl0dzWEa9s2RsvbtoT1RWl0d6Wys1b1Bz1TR1RVVkWxy+eFb936tJD3c55AgQIECBAgAABAgQIECBAgACBCS4gwDvBJ0j3CBypQLIq97JzT4yfPrUhNm2vS6VjaI7Wjq4or+iLpfNrUit3a+PNpy2LJF+vjQABAgQIECBAgAABAgQIECBAILcFBHhze/70nsCQAjVTK1JB3jWxc19TbNi6O7Zt3xnLly6OJfNmSMswpJiDBAgQIECAAAECBAgQIECAAIHcFBDgzc1502sChyUwt7Y6qsqLYnpZZ6xYMSfKyqzaPSw4lQgQIECAAAECBAgQIECAAAECOSLgJWs5MlG6SYAAAQIECBAgQIAAAQIECBAgQIAAgQMFBHgPFLFPgAABAgQIECBAgAABAgQIECBAgACBHBEQ4M2RidJNAgQIECBAgAABAgQIECBAgAABAgQIHCggwHugiH0CBAgQIECAAAECBAgQIECAAAECBAjkiIAAb45MlG4SIECAAAECBAgQIECAAAECBAgQIEDgQAEB3gNF7BMgQIAAAQIECBAgQIAAAQIECBAgQCBHBAR4c2SidJMAAQIECBAgQIAAAQIECBAgQIAAAQIHCgjwHihinwABAgQIECBAgAABAgQIECBAgAABAjkiIMCbIxOlmwQIECBAgAABAgQIECBAgAABAgQIEDhQQID3QBH7BAgQIECAAAECBAgQIECAAAECBAgQyBEBAd4cmSjdJECAAAECBAgQIECAAAECBAgQIECAwIECArwHitgnQIAAAQIECBAgQIAAAQIECBAgQIBAjggI8ObIROkmAQIECBAgQIAAAQIECBAgQIAAAQIEDhQQ4D1QxD4BAgQIECBAgAABAgQIECBAgAABAgRyRECAN0cmSjcJECBAgAABAgQIECBAgAABAgQIECBwoIAA74Ei9gkQIECAAAECBAgQIECAAAECBAgQIJAjAgK8OTJRukmAAAECBAgQIECAAAECBAgQIECAAIEDBQR4DxSxT4AAAQIECBAgQIAAAQIECBAgQIAAgRwREODNkYnSTQIECBAgQIAAAQIECBAgQIAAAQIECBwoIMB7oIh9AgQIECBAgAABAgQIECBAgAABAgQI5IiAAG+OTJRuEiBAgAABAgQIECBAgAABAgQIECBA4EABAd4DRewTIECAAAECBAgQIECAAAECBAgQIEAgRwQEeHNkonSTAAECBAgQIECAAAECBAgQIECAAAECBwoI8B4oYp8AAQIECBAgQIAAAQIECBAgQIAAAQI5IiDAmyMTpZsECBAgQIAAAQIECBAgQIAAAQIECBA4UKDkwAP2CRytQEdHR3R1dR3t5a7LkkBbW1v6zq2trdHZ2ZmlVtyWwNgItLe3p2+UPK/FxcVjc1N3IZAlgf7ntaWlJSZN8jPzLDG77RgJJH9PS7bkeS0qKhqju7oNgewI9D+vzc3N2WnAXQmMoYDndQwx3SrrAklMoK+vL/z5mnXqo2qgsrLyqP9dIcB7VOQuGkqgp6dHAHEomHE+1t3dne5B8gd5b2/vOPdG8wRGFhj4vArwjmzl7PgLDHxeBXjHfz70YGSBgc+rAO/IVs6Ov0Dy74pkszhh/OdCDw4tkDyvScDM83poKzXGX8DzOv5zMFIPKioqRjo94jkB3hF5nDwSgeQnDcnHNrEEkp/M1dfXx/Tp06OsrGxidU5vCBwg0NjYGA0NDenntbS09ICzdglMLIHkWU2e2ZqaGivOJ9bU6M0QAklQt6mpKf28+oHEEEAOTTiB5HmdMWPGhOuXDhE4UCAJ7ia/HeF5PVDG/kQUSAK8yW+heV4n4uyMrk9+n3B0fq4mQIAAAQIECBAgQIAAAQIECBAgQIDAuAkI8I4bvYYJECBAgAABAgQIECBAgAABAgQIECAwOgEB3tH5uZoAAQIECBAgQIAAAQIECBAgQIAAAQLjJiDAO270GiZAgAABAgQIECBAgAABAgQIECBAgMDoBAR4R+fnagIECBAgQIAAAQIECBAgQIAAAQIECIybgADvuNFrmAABAgQIECBAgAABAgQIECBAgAABAqMTEOAdnZ+rCRAgQIAAAQIECBAgQIAAAQIECBAgMG4CArzjRq9hAgQIECBAgAABAgQIECBAgAABAgQIjE5AgHd0fq4mQIAAAQIECBAgQIAAAQIECBAgQIDAuAkI8I4bvYYJECBAgAABAgQIECBAgAABAgQIECAwOgEB3tH5uZoAAQIECBAgQIAAAQIECBAgQIAAAQLjJiDAO270GiZAgAABAgQIECBAgAABAgQIECBAgMDoBAR4R+fnagIECBAgQIAAAQIECBAgQIAAAQIECIybgADvuNFrmAABAgQIECBAgAABAgQIECBAgAABAqMTEOAdnZ+rCRAgQIAAAQIECBAgQIAAAQIECBAgMG4CArzjRq9hAgQIECBAgAABAgQIECBAgAABAgQIjE5AgHd0fq4mQIAAAQIECBAgQIAAAQIECBAgQIDAuAkI8I4bvYYJECBAgAABAgQIECBAgAABAgQIECAwOgEB3tH5uZoAAQIECBAgQIAAAQIECBAgQIAAAQLjJiDAO270GiZAgAABAgQIECBAgAABAgQIECBAgMDoBAR4R+fnagIECBAgQIAAAQIECBAgQIAAAQIECIybgADvuNFrmAABAgQIECBAgAABAgQIECBAgAABAqMTEOAdnZ+rCRAgQIAAAQIECBAgQIAAAQIECBAgMG4CArzjRq9hAgQIECBAgAABAgQIECBAgAABAgQIjE5AgHd0fq4mQIAAAQIECBAgQIAAAQIECBAgQIDAuAkI8I4bvYYJECBAgAABAgQIECBAgAABAgQIECAwOgEB3tH5uZoAAQIECBAgQIAAAQIECBAgQIAAAQLjJiDAO270GiZAgAABAgQIECBAgAABAgQIECBAgMDoBAR4R+fnagIECBAgQIAAAQIECBAgQIAAAQIECIybgADvuNFrmAABAgQIECBAgAABAgQIECBAgAABAqMTEOAdnZ+rCRAgQIAAAQIECBAgQIAAAQIECBAgMG4CArzjRq9hAgQIECBAgAABAgQIECBAgAABAgQIjE5AgHd0fq4mQIAAAQIECBAgQIAAAQIECBAgQIDAuAkI8I4bvYYJECBAgAABAgQIECBAgAABAgQIECAwOgEB3tH5uZoAAQIECBAgQIAAAQIECBAgQIAAAQLjJiDAO270GiZAgAABAgQIECBAgAABAgQIECBAgMDoBAR4R+fnagIECBAgQIAAAQIECBAgQIAAAQIECIybgADvuNFrmAABAgQIECBAgAABAgQIECBAgAABAqMTEOAdnZ+rCRAgQIAAAQIECBAgQIAAAQIECBAgMG4CArzjRq9hAgQIECBAgAABAgQIECBAgAABAgQIjE5AgHd0fq4mQIAAAQIECBAgQIAAAQIECBAgQIDAuAkI8I4bvYYJECBAgAABAgQIECBAgAABAgQIECAwOgEB3tH5uZoAAQIECBAgQIAAAQIECBAgQIAAAQLjJiDAO270GiZAgAABAgQIECBAgAABAgQIECBAgMDoBAR4R+fnagIECBAgQIAAAQIECBAgQIAAAQIECIybgADvuNFrmAABAgQIECBAgAABAgQIECBAgAABAqMTEOAdnZ+rCRAgQIAAAQIECBAgQIAAAQIECBAgMG4CArzjRq9hAgQIECBAgAABAgQIECBAgAABAgQIjE5AgHd0fq4mQIAAAQIECBAgQIAAAQIECBAgQIDAuAkI8I4bvYYJECBAgAABAgQIECBAgAABAgQIECAwOgEB3tH5uZoAAQIECBAgQIAAAQIECBAgQIAAAQLjJlAybi1PwIZ7enriRz/6Uaxbty62bNkSDQ0NUVtbG8uXL4+LL744li1blrVeH8u2v/3tb8fPfvaz9FiuuuqqOO2007I2LjcmQIAAAQIECBAgQIAAAQIECBAgQCB7AgK8v7N9/PHH44tf/GJs3LjxIO21a9fGXXfdFe9+97vj2muvPej8aA8cy7aT4HUyziSgnGyXXnrpaLvvegIECBAgQIAAAQIECBAgQIAAAQIExklAgDcFnwRYP/zhD2eCnslcFBUVxeTJk6OtrS0zNXfffXc0NTXFjTfeGJMmjU12i2PZdjKWT37yk4PGmRmcAgECBAgQIECAAAECBAgQIECAAAECOSdQ8AHeDRs2xMc+9rFM0HP27Nlxww03xCmnnBKVlZXx0ksvxb333hsPPvhgenK/973vpevefPPNo57sY912snJ369ato+63GxAgQIAAAQIECBAgQIAAAQIECBAgMDEECj7Ae8cdd2RW6S5atChuvfXWSIK8/dvq1asj+SS5eO+888704STYe/XVV8f8+fP7qx3V97Fs+9FHH4377rvvqPrpIgIECBAgQIAAAQIECBAgQIAAAQIEJqbA2OQZmJhjO2SvNm/eHD/5yU8y9ZKVuwODu5kTqcIHPvCBOPvss9OH+vr64rvf/e7A00dcPpZt19XVxf/4H/8j3cck7URpaekR99cFBAgQIECAAAECBAgQIECAAAECBAhMPIGCDvA+8MADkQRrk23ZsmVxxhlnjDhDV155Zeb8/fffHx0dHZn9Iy0cy7b/9m//Nurr69Nd/LM/+7N0buEj7a/6BAgQIECAAAECBAgQIECAAAECBAhMPIGCDvA+++yzmRk566yzMuXhCkle3rKysvTpxsbGSNIeHO12rNr+zne+E4899li6m294wxvi8ssvP9ouu44AAQIECBAgQIAAAQIECBAgQIAAgQkmULAB3u7u7njhhRcy03HyySdnysMVktQGJ5xwQuZ08gK2o9mOVdsbN26M22+/Pd3FadOmxU033XQ03XUNAQIECBAgQIAAAQIECBAgQIAAAQITVKBgA7xJ8LOzszMzLYf7wrS5c+dmrknucTTbsWg7CSJ/8pOfzIzxL//yL9Mvijua/rqGAAECBAgQIECAAAECBAgQ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IECAAAECE0rgmAV4kxWwX/va19IrW/fu3Zt+oVlXV9dRYTz++ONHdd1wF1VXV8dnPvOZ+Od//ue47777Yt++fQdVPeWUU+JP//RP48QTTzzo3GgOjGfbo+m3awkQIECAAAECBAgQIECAAAECBAgQGH+Bor7Uls1u7Ny5Mz784Q/Ht771rXRQdyzaynKXIwlAv/jii+m0DPPnz49FixbFvHnzxqLrh7zHeLZ9yM6pkJMCzc3N6RQkK1askIM3J2ewsDqdvNBy27Ztcdxxx2VegllYAkabSwINDQ2xffv2OP7446O4uDiXuq6vBSiwf//+SP5efsIJJ8SkSQX9nuUCnP3cG3LyDpRdu3bF6tWrc6/zelxwAskCsWRB26pVqwpu7AacewLJs5r8nSD5+6stvwSyuoI3eWjOP//8dLA0l9hmzpwZyWc8tvFsezzGq00CBAgQIECAAAECBAgQIECAAAECBI5eIKsB3o9//OMjBneTHLRWEBz95LmSAAECBAgQIECAAAECBAgQIECAAIHCFshagLe7uzu+8pWvDNJNloAnuW7PPvvs9ArZkpKsNT+oXTsECBAgQIAAAQIECBAgQIAAAQIECBDIR4GsRVjXrVsX7e3tGbNzzz03/u3f/i1mzJiROaZAgAABAgQIECBAgAABAgQIECBAgAABAkcvkLU3LDz55JODevXlL39ZcHeQiB0CBAgQIECAAAECBAgQIECAAAECBAiMTiBrAd6nnnoq07N58+Z5o2RGQ4EAAQIECBAgQIAAAQIECBAgQIAAAQJjI5C1AO/A/LrnnHPO2PTWXQgQIECAAAECBAgQIECAAAECBAgQIEAgI5C1AO+CBQsyjWzfvj1TViBAgAABAgQIECBAgAABAgQIECBAgACBsRHIWoD3/PPPz/Twueeei56ensy+AgECBAgQIECAAAECBAgQIECAAAECBAiMXiBrAd7TTjstTjzxxHQPm5ub40tf+tLoe+sOBAgQIECAAAECBAgQIECAAAECBAgQIJARyFqAN2nhtttui6KionRjt9xyS7z66quZhhUIECBAgAABAgQIECBAgAABAgQIECBAYHQCWQ3wJmkaPv3pT6d72NDQECeddFJ8/OMfj7a2ttH12tUECBAgQIAAAQIECBAgQIAAAQIECBAgECXZMujq6ook9+6ll14au3fvjltvvTXa29vjE5/4RPzd3/1dLFq0KJYsWZL+VFdXH1E3knvZCBAgQIAAAQIECBAgQIAAAQIECBAgUOgCWQvwbt++PZI8vENtHR0d8corr6Q/Q50/1DEB3kMJOU+AAAECBAgQIECAAAECBAgQIECAQCEIZDVFQyEAGiMBAgQIECBAgAABAgQIECBAgAABAgTGS0CAd7zktUuAAAECBAgQIECAAAECBAgQIECAAIFRCmQtRcP8+fPjqaeeGmX3XE6AAAECBAgQIECAAAECBAgQIECAAAECwwlkLcBbWloar3vd64Zr13ECBAgQIECAAAECBAgQIECAAAECBAgQGKWAFA2jBHQ5AQIECBAgQIAAAQIECBAgQIAAAQIExktAgHe85LVLgAABAgQIECBAgAABAgQIECBAgACBUQqMe4C3t7c3tm7dGnV1ddHX1zfK4bicAAECBAgQIECAAAECBAgQIECAAAEChSOQtRy8wxE+88wzcc8998QLL7wQL774Yrz88svR1taWrl5SUhIzZ86MJUuWxOWXXx5XXHFFrFy5crhbOU6AAAECBAgQIECAAAECBAgQIECAAIGCFjhmAd7m5ub467/+6/jCF74Q3d3dQ6Inx3fu3Jn+PPbYY/HRj340Lrroovj6178e8+fPH/IaBwkQIECAAAECBAgQIECAAAECBAgQIFCoAsckRcN9990Xq1atis997nPDBneHm4CHHnooTjnllPjBD34wXBXHCRAgQCAPBLq6e2LnvubYsa8l9uxvid5eaXvyYFoNgQABAgQIECBAgAABAgSyLJD1FbzPPfdcvPOd74yOjo6DhlJcXBwLFy6MpUuXRmtra2zcuDH27NlzUL19+/bFO97xjvjlL38Zq1evPui8AwQIECCQuwKdXT3xq3Vb4sVNe2JPXWM0NjbF05taY2bNlDj5uHnpz6RJRbk7QD0nQIAAAQIECBAgQIAAAQJZFMhqgDdJuXD11VcPCu6Wl5fHddddF9dee206sFtaWjpoeEk+3iQ3b7La9+67706t4OpNn29paUkHipMgb0VFxaBr7BAgQIBAbgq0tnfG93/2Qry0aW/sbWiJ8pJJkQR899S3xNbdjbG7rjl2pT6/f9ZxUTzpmPzSSW5C6jUBAgQIECBAgAABAgQIFKxAVv+1/OlPfzqefPLJDO6ll14aL730Unz2s59NvzztwOBuUjEJ3p56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dyx9fnC608EV7BE4BgIzUyvNqypTKURSP5iYMkJu3eQHEXNnVqdSOVQdg15pggABAgQIECBAgAABAgQI5JZAVgO8W7dujYsvvjh2796dWyp6SyBHBaZXV6Rfrtbc0hFTJg//Mqrm1s50QO1QK31zlEG3c0TguEW1Mbe2Op5LvfRvWtXkIYO8e1P5pJNt0expqboCvDkytbpJgAABAgQIECBAgAABAsdQYOgknWPUgU9+8pNjFtx905veFF/60pfilVdeGaPeuQ2B/BNYvmBG+oVVSR7eru6eIQfY1t6VSs/Qka63aM70Ies4SOBYCFSmfghxzslLYum8mti4Y3/s3NcUnV090ZtK1ZOkGklewra3oSVWLKyNN526bNgXsR2LvmqDAAECBAgQIECAAAECBAhMVIGsreBtaWmJr3/964PGfdVVV8Wf/MmfxIoVK6K6ujqWL18e+/bti5KSkli/fn1UVFTE/v370/l2k5QO//RP/5QJELe3t8f73ve+KCsbflXioMbsEChAgWQF7ykr50VTagXvK1v2xaI506Ks5LVcDPVNbbFtd2Msnjs9zli9MJIcqDYC4ymwetns6Enl2P3FM5vSAd5NO5pSL1Nri+lTI5VLujId3D3/jBUxe4bVu+M5T9omQIAAAQIECBAgQIAAgYkrkLUA79q1a9MvUusfepKH9+///u/7d9Pf5513XnznO99JvzztscceiyuuuCJmzpwZK1eujPPPPz8+9KEPRRIU/sEPfhBPPPFE3HTTTfH5z39+0D3sECAwWOCsExenVj92xzMv74jtexqjvaMzuro6Y0d9dzp9w/KFM+LMNYvixBVzB19oj8A4CZyUehaXzK2JFzbtjk3b9sSevXWxaOG81A8iZsSqpbP9IGKc5kWzBAgQIECAAAECBAgQIJAbAllL0fDkk09mBKZNmxYf/ehHM/v9hQsuuKC/GA888ECm3F+YNWtWfPvb3441a9akD33xi1+Ml156qf+0bwIEhhCYNKkoLkiteHzrm1fHG1+3NFYunhlzZ1TEmuWz482nLYvLzzsxzjxx0RBXOkRg/ASqp5Snf/Bw8dnHxfmnLohLzjk+vRrdKvPxmxMtEyBAgAABAgQIECBAgEBuCGRtBW9dXV1G4JxzzomamprMfn/hda97XX8xfvWrX2XKAwtTpkxJ595NVvT29PTEpz/96fjf//t/D6yiTIDAEALLU3lLk09DQ2Os37gpVp2wMiomTx6ipkMECBAgQIAAAQIECBAgQIAAAQK5KpC1Fbz19fUZk2XLlmXKAwsnnHBCZveFF14YlNIhcyJVSFI5HH/88elDDz300MBTygQIHEKguHhSTC4rieJJWfu/+yF64DQBAgQIECBAgAABAgQIECBAgEC2BLIW8Uly6fZvSaqFobbZs2dHkr4h2bq7u2PdunVDVUsfe/Ob35z+3rFjR/rFbMNWdIIAAQIECBAgQIAAAQIECBAgQIAAAQIFIpC1AO+qVasyhHv27MmUDywMXMX79NNPH3g6s7948eJM+ZlnnsmUFQgQIECAAAECBAgQIECAAAECBAgQIFCoAsckwLt58+ZhfVeuXJk5N1KAt7e3N1NPgDdDoUCAAAECBAgQIECAAAECBAgQIECAQAELZC3AO3Bl7g9/+MPYuXPnkMwDA7w///nPh6yTHFy/fn3mXHFxcaasQIAAAQIECBAgQIAAAQIECBAgQIAAgUIVyFqAd+rUqTF//vy0a2dnZ1x//fVDvkTttNNOy9j/4he/iOeffz6z319IXth2//339+/GihUrMmUFAgQIECBAgAABAgQIECBAgAABAgQIFKpA1gK8Ceg111yTcb3nnnvivPPOi+9+97vpF6r1n7jgggsiCQYnW5KG4fLLL4+6urr+09Ha2hrvf//7Y//+/Zljxx13XKasQIAAAQIECBAgQIAAAQIECBAgQIAAgUIVyGqA94Ybbohp06ZlbJMVum9/+9vj2muvzRyrrq6O//Jf/ktm/+WXX46FCxfGW9/61njPe94Tycva/s//+T+Z8yeddFIsW7Yss69AgMDQAs1tnfHE81vj/z72Sjzy1Lb4v794OX7z0vbo6Owe+gJHCRAgQIAAAQIECBAgQIAAAQIEck6gJJs9rq2tjXvvvTfe8pa3REdHR6apJIA7cPvQhz4Ut912W7S1taUPJ98PPPDAwCqZ8j/8wz/EpElZjUtn2lIgkKsCL2zcHWt/szF27muKuoaWaG5pjY2726Nm2q54+uUdcf4ZK2LRnOm5Ojz9JkCAAAECBAgQIECAAAECBAgQ+J1A1iOl559/fjz88MNx6qmnZtAPzKG7YMGCSF7ENmPGjEydoQrvfe9745JLLhnqlGMECPxO4KXNe+JHj78Sz2/YnTpSFIvnTIulc6pjYeq7s6snnn11Z/x/a1+MHXsamREgQIAAAQIECBAgQIAAAQIECOS4QNYDvInPG97whnjiiSfiwQcfjA984AORpFk4cHvjG98Yjz76aJxzzjlRXFw86HRVVVV87WtfizvuuGPQcTsECAwWaE+lX/j505ti4/a6WDR3esybWR3lZSWpVe9FUVFeml61O6umKl7dui9+lqrX29s3+Ab2CBAgQIAAAQIECBAgQIAAAQIEckogqykaBkokaRUuvfTS9Gfg8YHl1atXx89//vNobGxMB3t37NiRXvmbBITLy8sHVlUmQGAIgfWpwO3OvU0xdcrkqK4c+v8ztdMqo76pLbbuqo9texqkahjC0SECBAgQIECAAAECBAgQIECAQK4IHLMA75GATJ06Nf79v//3R3KJugQIpAT27G+O5taOmDWjakSPaVWT0/X21rcI8I4o5SQBAgQIECBAgAABAgQIECBAYGILHJMUDRObQO8I5I9AV09v9KTSLhQf4kWEJcWTUvV6oyuVk9dGgAABAgQIECBAgAABAgQIECCQuwICvLk7d3pO4CCBqoqyKC8tjo5ULt6RtraOrlS9kpgyTBqHka51jgABAgQIECBAgAABAgQIECBAYOIIHHWKhnvvvTceeuihzEiWL18eN954Y2a/rq4ubrnllsz+WBZuv/32sbydexHIG4GFs6dHzdTK2L63MaZXV0RR0cFDS1b47m9si+MWzYyFs6cdXMERAgQIECBAgAABAgQIECBAgACBnBE46gDv2rVr40tf+lJmoG984xsHBXibmpoGnc9UHIOCAO8YILpFXgosnDMtli2YEXWNrbF55/6D8ut2p1I4bNy+PxUErohVS2dHkovXRoAAAQIECBAgQIAAAQIECBAgkLsCRx3gzd0h6zmB/BY4/4wV0dLWGS9v3hPPb9wdU8pLoqOjPVq7G6KptTNmTZ8Sq5fNiTecsiS/IYyOAAECBAgQIECAAAECBAgQIFAAAgK8BTDJhlhYAtWpvLqXnbcmHntmc7y0eW/UNTRHKr4blZNLY9HcmnRw98w1C6O0pLiwYIyWAAECBAgQIECAAAECBAgQIJCHAkcd4E3y7V5zzTUZksrKykw5KcyfPz+ee+65QcfsECBwbAQqJ5fFBWceF2edtDg2bNkVm7duj5UrlqYCvDNicnnpsemEVggQIECAAAECBAgQIECAAAECBLIucNQB3jlz5kTyGW4rLS2NNWvWDHfacQIEjoHAlIqyWDq/Jop7mmPJvJooKxPcPQbsmiBAgAABAgQIECBAgAABAgQIHDOBScesJQ0RIECAAAECBAgQIECAAAECBAgQIECAwJgKZC3Au2PHjrj22mvjscceG9MOuxkBAgQIECBAgAABAgQIECBAgAABAgQI/FYgawHezs7O+PKXvxznnHNOOlXDZz7zmUiCvjYCBAgQIECAAAECBAgQIECAAAECBAgQGBuBrAV4B3bv+eefj5tuuikWLVoUb33rW+Pb3/52JAFgGwECBAgQIECAAAECBAgQIECAAAECBAgcvcAxCfD2d6+npyceeOCBuOKKK2L+/PnxoQ99KJ566qn+074JECBAgAABAgQIECBAgAABAgQIECBA4AgEshbgXbx4cTzyyCPxvve9L6ZOnXpQl/bt2xdf+MIX4rTTTkt/kvLevXsPqucAAQIECBAgQIAAAQIECBAgQIAAAQIECAwtkLUAb1FRUZx33nnx1a9+NXbt2hXf/OY34w/+4A+ipKTkoJ4kq3iT1bwLFixIr+69//77I1ntayNAgAABAgQIECBAgAABAgQIECBAgACB4QWyFuAd2OTkyZPjXe96VySB223btsX//J//M04//fSBVdLlJC9vkp/3bW97WyxcuDA+/OEPR5K/10aAAAECBAgQIECAAAECBAgQIECAAAECBwscvJz24DpjemT27Nlx/fXXpz9J8PbOO++Mu+66K7Zs2TKonZ07d8bf//3fpz9nn312/PEf/3H84R/+YUybNm1QvbHcSVYN/+hHP4p169al+9PQ0BC1tbWxfPnyuPjii2PZsmVj2dyge2Wj7WTl9P/7f/8vNm3aFJs3b46+vr5I/E8++eT4/d///aipqRnUBzsECBAgQIAAAQIECBAgQIAAAQIECOSWQFEq6Nc33l1OupDk602CvckK3sbGxiG7VFFREf/hP/yHuPvuu4c8P5qDjz/+eHzxi1+MjRs3Dnubd7/73XHttdcOe/5oT4x1221tbemgeZIWI1kVPdRWXFycDpgnOZJLS0uHquJYngg0Nzenf2CxYsWKKCsry5NRGUa+CiR//ie/6XHcccf5sylfJzkPxrWvviWe37g7Nm/fG3v37Y8li+bHork1sXrZ7Kgo99/UPJjivBzC/v37I1lA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CCCAAAIIIIAAAgjME4GQBHj37Nmjd77znaqtrR2TLcasJm1TLhQUFDiPnJwc9fT0yObibWtr04EDB+R2u497bVdXl6644grZPIw2dQMFAQQQQGBuCsTERGuJycmbmx6nWhPTXbasRORln5vv9VwYVUx0tMoKMpWRFKW61CGtWFEk+7sOBQEEEEAAAQQQQAABBBCYCYEpB3irqqq0detWtbe3j+i//cP8oosu0tVXX60LLrhg3D/UBwYGtHPnTj355JO69957VVdX57+Wx+PR3//93ysqKkrM5PWzsIEAAggggAACCCCAAAIIIIAAAggggAACCGjKOXivvfba44K773rXu1RTU6NHH31UF1988bjBXfseJCQk6Oyzz9Ytt9yiiooKffe731VaWpr/7fEFeaurq/3H2EAAAQQQQAABBBBAAAEEEEAAAQQQQAABBOa7wJQCvL/85S/129/+doThjTfeqMcff1z5+fkjjk90x878tQu0vfLKK1q6dKn/Zb4F3PwH2EAAAQQQQAABBBBAAAEEEEAAAQQQQAABBOa5wJQCvPfff/8Ivve///268847FW1y0021rFixQvfcc8+Iy9iZvTadAwUBBBBAAAEEEEAAAQQQQAABBBBAAAEEEEBAk0/RYBdJe+KJJ/yG2dnZ+spXvuLfD8XG+eefr0suucR/qaamJj3yyCP+fTYQQAABBBBAAAEEEEAAAQQQQAABBBBAAIH5LDDpqbZ79+5VX1+f3+6KK65QXl6efz9UG7fffvuIS9nUDRQEEEAAAQQQQAABBBBAAAEEEEAAAQQQQACBKczgPXLkyAi/s846a8R+qHbWrVunzMxM/+Vqa2v922wggAACCCCAAAIIIIAAAggggAACCCCAAALzWWDSM3gbGhpGuK1evXrEfqh2oqKitH79ev/lCPD6KdhAAAEEEEAAAQQQQAABBBBAAAEEEEAAgXkuMOkAb2B6BmuYlpYWNsqsrCz/tQnw+inYQAABBBBAAAEEEEAAAQQQQAABBBBAAIF5LhARAd7A4HFdXd08f8sYPgIIIIAAAggggAACJycwNDSsQe/Qyb2I2ggggAACCCCAAAIRIRA72V663e4RL01KShqxH8qdwGuPbjeU7XAtBBBAAAEEEEAAAQTmisDw8LCO1rTqQEWTKmtb1OlyaVdFn8oKs7R+eaGy05PnylAZBwII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66+/rg0bNvhfF7hRVlbm333jjTechdj8B0K0YVNA/OEPf3ACzdXV1ebWy07l5ORoyZIlevvb365gM5FD0Xw42m5padGf//xn1dbWyo7HBsaLiopUWlqqRYsWyQbNExMTQ9F9roEAAggggAACCCCAAAIIIIAAAgggMEsE9lc0qdoEcBPiYlWYm35cr1KS4rWoMEtVps5rB2q13AR9o6KijqvHgcgQmJYAb1VVVVCN5cuXjwjwBqs4ZG7d9ZVwBHhtkNku4lZRUeFrxv+8Y8cO/fjHP9aVV16pj33sY/7jodoIddv9Jp/KT37yEz344IOy24Hl4MGD/l07g/r666/X1q1b/cfYQAABBBBAAAEEEJhbAu2uXn3kzsed9Sw2rKrWrde+fW4NkNEggAACCCCAwHECNrjbZn4HWFSUfdw53wEb5I2LjXZSNbR29io3M8V3iucIEwhbiobA/LpPPvmkGhoaxqSxAV5fef75532bxz0fO3bMfyzUC5nYAOunP/3pEcFd+61FUlKSv0278cADD+i///u/FRhsHlFhEjuhbru3t1fXXHON7rvvvhHBXTuetLS0ET1samrSzTffrC984QsjjrODAAIIIIAAAgggMHcEvEPDau7oVVvXgNpdfXNnYIwEAQQQQAABBIIKdPX0y6ZpSIwff25nUkKcBky97l530GtxYvYLjP8uT6H/6enpTjoAu7q02+3Wdddd58yCHb2I2mmnneZv5YUXXtD+/fu1evVq/zG7YRdse+KJJ/zHli5d6t+e6kZ5ebk++9nPOjMa7LXsrNZPfvKTWr9+vZKTk52F3R577DH9+te/dpp6/PE3Zz985jOfmWrTCkfbX/7ylxU4Y/rss8/WBz/4QSe9REJCglm8yKXXXntNd999tz/obse2ceNGveMd75jymLgAAggggAACCCCAAAIIIIAAAggggMDMCsTGRDspF4bMF73R0cFTLzjnzWlbnxK5AmF9966++mq/zMMPP6xt27bpF7/4hbOgmu/EueeeKxsMtsXOjL3kkkvU1tbmOy07I/UjH/mI2tvb/ceWLVvm357qxve+9z319b05k8Hmpv3Wt76lc845x5ntamcK22DzjTfeqKuuusrflA2I2sD1VEuo23766af1u9/9zt+tD3/4w/riF78omw/ZBndtsdb2fbj//vsVOMv6K1/5imzOXgoCCCCA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CIQT6dFHH/X7QP7FL34hF198sXTt2lUaNWokhYWFxuOa+iGvj79pmjt3rhx22GHCI5lOXCF31+H0/fr444/LqlWrLBR97PKSSy4xhhLRHyX0l239seLdd981ynz55Zei9/htt91m7UPGuwJO3691SesTGj/99FNdRdjmcYFY3K/fffedpa7DMzVr1sxaJoNAXQJO36/64+3vf/97Y0gxbZf2zNVhbs4880xjaDMdosH8DmB+p549e7bxOLF+dyB5W8Cp+1U7H/Tr1y8sbO188+OPP1r7HHzwwX6dGKwNZDwj4NT9aoIuWbJEXnvtNXPRGO7muuuuMz4/tSf57t275R//+IfR4UY7KuiPaPfdd588++yzPCVhqTmbIcDrrDe1IRCSgAalNIClY9s4kebPny//+c9/rKp+97vfyUUXXWQta0b/uDv22GNlyJAhctVVVxnjmun6P//5z3L44YfzOIZieDQ5fb+uXr1aZs2aZWlfeOGFctlll1nLmtExo++44w7jvtQevJr+9a9/Gb8s6zh8JO8KOH2/1iW9cOFCefvtt+sqwjaPC8TiftUhb8xAmPLffffd0rt3b49fCU4/FIFY3K9///vf/YK7+kTEL3/5S6u5GoQ45JBDjOEb9PFh85Hip59+WnRIEf1RmORNASfv1x49ehhBr1CldT6JW2+91QrwtmrVSh544AG/HumhHotyiSHg5P1qij388MNm1hjHXIfB8Z2HR3vxageb7t27G393aWEd1vGNN94whsizdibjmABDNDhGTUUI1C+gk0XddNNNMmnSJMeCu9oqDfCaSf+I08crgiX9UL/rrrusLxj663JDJmIJdnzWx49ArO5XHQfKTDqwv/Y0D5Z07F3fR4j0yxHJmwKxul+DaWsvjIceesjYrJ+rjKkXTMqb62N5v2qPMXPIEO0N2a1bN29eBM46ZIFY3a/6qLLOG2GmK664wi+4a67Xd30aTeePMJOOJclQJKaGt95jdb+Go6z3tTnxqn4O69+GjG0ajmDilI3V/ZqXlydat5m0565vcNdcr+/aCWzo0KHWKmIDFoXjGQK8jpNTIQKBBXQmTH2k7NNPP7UKtG/f3njEzFphU8b3Q/jUU0+1grfBqtMekDp0g5l8J2Az1/Ge2AKxul91EhUd48lMxx13XJ2PAGnQzLfn2YYNG8xdefeQQKzu17qIH3zw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" />

<!-- rnb-plot-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



## Calculate the performance of each pipeline

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzIDwtIGFjY3VyYWN5X3Njb3JlcyAlPiUgZmlsdGVyKHRlc3Rpbmdfc2V0X2luZGV4PD01KVxuYGBgIn0= -->

```r
accuracy_scores <- accuracy_scores %>% filter(testing_set_index<=5)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzJG51bXJ1bGVzIDwtIHN0cl9jb3VudChhY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvLCBcIjFcIilcbiNhY2N1cmFjeV9zY29yZXMgPC0gYWNjdXJhY3lfc2NvcmVzW29yZGVyKGFjY3VyYWN5X3Njb3JlcyRudW1ydWxlcywgZGVjcmVhc2luZz1GKSxdXG5hY2N1cmFjeV9zY29yZXMgPC0gYWNjdXJhY3lfc2NvcmVzW29yZGVyKGFjY3VyYWN5X3Njb3JlcyRNQ0MsIGRlY3JlYXNpbmc9RiksXVxuYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibyA8LSBmYWN0b3IoYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibywgbGV2ZWxzID0gdW5pcXVlKGFjY3VyYWN5X3Njb3JlcyR0b29sY29tYm8pKVxuYGBgIn0= -->

```r
accuracy_scores$numrules <- str_count(accuracy_scores$toolcombo, "1")
#accuracy_scores <- accuracy_scores[order(accuracy_scores$numrules, decreasing=F),]
accuracy_scores <- accuracy_scores[order(accuracy_scores$MCC, decreasing=F),]
accuracy_scores$toolcombo <- factor(accuracy_scores$toolcombo, levels = unique(accuracy_scores$toolcombo))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


accounting for rules that have multiple tools in them

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzJG51bXRvb2xzIDwtIGFjY3VyYWN5X3Njb3JlcyRudW1ydWxlc1xuXG5mb3IgKGkgaW4gMTpucm93KGFjY3VyYWN5X3Njb3JlcykpIHtcbiAgdG9vbGNvbWJvIDwtIGFzLmNoYXJhY3RlcihhY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvW2ldKVxuICBpZiAoc3Vic3RyKHRvb2xjb21ibywgMSwgMSk9PVwiMVwiKSB7XG4gICAgYWNjdXJhY3lfc2NvcmVzJG51bXRvb2xzW2ldIDwtIGFjY3VyYWN5X3Njb3JlcyRudW10b29sc1tpXSArIDNcbiAgfVxuICBcbiAgaWYgKHN1YnN0cih0b29sY29tYm8sIDUsIDUpPT1cIjFcIikge1xuICAgIGFjY3VyYWN5X3Njb3JlcyRudW10b29sc1tpXSA8LSBhY2N1cmFjeV9zY29yZXMkbnVtdG9vbHNbaV0gKyA0XG4gIH1cbn1cblxuYWNjdXJhY3lfc2NvcmVzJG51bXRvb2xzW2FjY3VyYWN5X3Njb3JlcyRudW10b29scz42XSA8LSA2XG5gYGAifQ== -->

```r
accuracy_scores$numtools <- accuracy_scores$numrules

for (i in 1:nrow(accuracy_scores)) {
  toolcombo <- as.character(accuracy_scores$toolcombo[i])
  if (substr(toolcombo, 1, 1)=="1") {
    accuracy_scores$numtools[i] <- accuracy_scores$numtools[i] + 3
  }
  
  if (substr(toolcombo, 5, 5)=="1") {
    accuracy_scores$numtools[i] <- accuracy_scores$numtools[i] + 4
  }
}

accuracy_scores$numtools[accuracy_scores$numtools>6] <- 6
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzJG51bXJ1bGVzIDwtIGFzLmZhY3RvcihhY2N1cmFjeV9zY29yZXMkbnVtcnVsZXMpXG5hY2N1cmFjeV9zY29yZXMkbnVtdG9vbHMgPC0gYXMuZmFjdG9yKGFjY3VyYWN5X3Njb3JlcyRudW10b29scylcbmBgYCJ9 -->

```r
accuracy_scores$numrules <- as.factor(accuracy_scores$numrules)
accuracy_scores$numtools <- as.factor(accuracy_scores$numtools)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzJHJ1bGV0eXBlIDwtIFwib3RoZXJcIlxuYWNjdXJhY3lfc2NvcmVzJHJ1bGV0eXBlW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm8gJWluJSBjKFwiMCAwIDEgMCAwIDFcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcIjAgMCAxIDEgMCAwXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCIwIDAgMCAxIDAgMFwiKV0gPC0gXCJoaWdoIHByZWNpc2lvblwiXG5hY2N1cmFjeV9zY29yZXMkcnVsZXR5cGVbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibyAlaW4lIGMoXCIxIDEgMCAwIDEgMVwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiMSAxIDAgMSAxIDFcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcIjEgMSAwIDAgMCAxXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCIxIDEgMCAxIDAgMVwiKV0gPC0gXCJoaWdoIHJlY2FsbFwiXG5hY2N1cmFjeV9zY29yZXMkcnVsZXR5cGVbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibyAlaW4lIGMoXCIwIDAgMSAwIDAgMVwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiMCAxIDEgMCAwIDFcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcIjAgMCAwIDAgMCAxXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCIxIDAgMSAxIDAgMVwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiMSAwIDEgMCAxIDFcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcIjEgMCAxIDAgMCAxXCIpXSA8LSBcImhpZ2ggTUNDXCJcbiNhY2N1cmFjeV9zY29yZXMkcnVsZXR5cGVbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibyAlaW4lIGMoXCIwIDAgMSAwIDAgMVwiKV0gPC0gXCJoaWdoIE1DQyBhbmQgaGlnaCBwcmVjaXNpb25cIlxuYWNjdXJhY3lfc2NvcmVzJHJ1bGV0eXBlW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm8gJWluJSBjKFwiMSAxIDEgMSAxIDFcIildIDwtIFwiYWxsXCJcblxuXG5gYGAifQ== -->

```r
accuracy_scores$ruletype <- "other"
accuracy_scores$ruletype[accuracy_scores$toolcombo %in% c("0 0 1 0 0 1",
                                                          "0 0 1 1 0 0",
                                                          "0 0 0 1 0 0")] <- "high precision"
accuracy_scores$ruletype[accuracy_scores$toolcombo %in% c("1 1 0 0 1 1",
                                                          "1 1 0 1 1 1",
                                                          "1 1 0 0 0 1",
                                                          "1 1 0 1 0 1")] <- "high recall"
accuracy_scores$ruletype[accuracy_scores$toolcombo %in% c("0 0 1 0 0 1",
                                                          "0 1 1 0 0 1",
                                                          "0 0 0 0 0 1",
                                                          "1 0 1 1 0 1",
                                                          "1 0 1 0 1 1",
                                                          "1 0 1 0 0 1")] <- "high MCC"
#accuracy_scores$ruletype[accuracy_scores$toolcombo %in% c("0 0 1 0 0 1")] <- "high MCC and high precision"
accuracy_scores$ruletype[accuracy_scores$toolcombo %in% c("1 1 1 1 1 1")] <- "all"

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzX21lbHQgPC0gYWNjdXJhY3lfc2NvcmVzICU+JSBcbiAgc2VsZWN0KHRlc3Rpbmdfc2V0X2luZGV4LCBwcmVjaXNpb24sIHJlY2FsbCwgTUNDLCBwcm9wX3ZpcmFsLCBudW1ydWxlcywgbnVtdG9vbHMsIHRvb2xjb21ibywgcnVsZXR5cGUpICU+JVxuICBwaXZvdF9sb25nZXIoY29scz1jKHByZWNpc2lvbiwgcmVjYWxsLCBNQ0MsIHByb3BfdmlyYWwpLCBcbiAgICAgICAgICAgICAgIG5hbWVzX3RvPVwicGVyZm9ybWFuY2VfbWV0cmljXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJwZXJmb3JtYW5jZV9tZXRyaWNfc2NvcmVcIilcbmBgYCJ9 -->

```r
accuracy_scores_melt <- accuracy_scores %>% 
  select(testing_set_index, precision, recall, MCC, prop_viral, numrules, numtools, toolcombo, ruletype) %>%
  pivot_longer(cols=c(precision, recall, MCC, prop_viral), 
               names_to="performance_metric",
               values_to="performance_metric_score")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZmlnIDwtIGdncGxvdChhY2N1cmFjeV9zY29yZXNfbWVsdFthY2N1cmFjeV9zY29yZXNfbWVsdCRwZXJmb3JtYW5jZV9tZXRyaWM9PVwicHJvcF92aXJhbFwiLF0sIFxuICAgICAgICAgICAgICBhZXMoeD1ydWxldHlwZSwgeT1wZXJmb3JtYW5jZV9tZXRyaWNfc2NvcmUsIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGNvbG9yPXJ1bGV0eXBlLCBmaWxsPXJ1bGV0eXBlKSkgK1xuICBnZW9tX2JveHBsb3QoKSArXG4gICNnZW9tX3BvaW50KCkgK1xuICB0aGVtZV9saWdodCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBheGlzLnRpY2tzLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIixcbiAgICBheGlzLnRleHQueT1lbGVtZW50X3RleHQoc2l6ZT0xMCksXG4gICAgYXhpcy50ZXh0Lng9ZWxlbWVudF90ZXh0KHNpemU9MTIsIGFuZ2xlID0gOTAsIHZqdXN0ID0gLTAuMSksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIHN0cmlwLmJhY2tncm91bmQgPSBlbGVtZW50X3JlY3QoZmlsbD1cIndoaXRlXCIsIGNvbG9yPVwiZ3JleVwiKSxcbiAgICBzdHJpcC50ZXh0ID0gZWxlbWVudF90ZXh0KGNvbG9yPVwiYmxhY2tcIiwgc2l6ZT0xMilcbiAgKSArXG4gIHhsYWIoXCJFbnZpcm9ubWVudFwiKSArXG4gIHlsYWIoXCJQcm9wb3J0aW9uIENhbGxlZCBWaXJhbFwiKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKGMocmV2KG1hZ21hKDcpWzM6Nl0pLCBcImdyZXlcIiksIDAuNCkpICtcbiAgc2NhbGVfY29sb3JfbWFudWFsKG5hbWU9XCJcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKGMocmV2KG1hZ21hKDcpWzM6Nl0pLCBcImdyZXlcIiksIDAuNikpICtcbiAgZ2VvbV9obGluZSh5aW50ZXJjZXB0ID0gMC4xLCBsaW5ldHlwZT1cImRhc2hlZFwiLCBjb2xvcj1cImRhcmsgZ3JleVwiKVxuXG5maWdcbmBgYCJ9 -->

```r
fig <- ggplot(accuracy_scores_melt[accuracy_scores_melt$performance_metric=="prop_viral",], 
              aes(x=ruletype, y=performance_metric_score, 
                                  color=ruletype, fill=ruletype)) +
  geom_boxplot() +
  #geom_point() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=10),
    axis.text.x=element_text(size=12, angle = 90, vjust = -0.1),
    legend.text=element_text(size=12),
    axis.title=element_text(size=12),
    strip.background = element_rect(fill="white", color="grey"),
    strip.text = element_text(color="black", size=12)
  ) +
  xlab("Environment") +
  ylab("Proportion Called Viral") +
  scale_fill_manual(name="",
                     values = alpha(c(rev(magma(7)[3:6]), "grey"), 0.4)) +
  scale_color_manual(name="",
                     values = alpha(c(rev(magma(7)[3:6]), "grey"), 0.6)) +
  geom_hline(yintercept = 0.1, linetype="dashed", color="dark grey")

fig
```

<!-- rnb-source-end -->

<!-- rnb-plot-begin eyJoZWlnaHQiOjQzMi42MzI5LCJ3aWR0aCI6NzAwLCJzaXplX2JlaGF2aW9yIjowLCJjb25kaXRpb25zIjpbXX0= -->

<img src="data:image/png;base64,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IBT915gzjB0B/b27qKjInPC9995bbXDXktHfyfv06SMrV6407+PWdp4RqGkBArw1fQXoHwEEEEAgpgVatmxpboamN0Tbvn27CfRqsFe/wq05/l599VXzaNu2rQn06uxeDQRREEAgcgX065yjRo2K3BNg5AjEqAABoRi98Jx2VAro791aNO+ur6nN9Bs5GuDV39UpCISLQNXzz8NllIwDAQQQQACBGBBo3769TJgwQd5++22ZO3euDB8+XM4991xz5jt37jQz/zQo9L//+7+mjs72pSCAAAIIIIBAaAQCDQjpKAkIheZa0QsC1QmcOXPGVElKSvJ69q7VZn5+vllMTk62NvGMQI0LMIO3xi8BA0AAAQQQQKCygN5sTR/33XefrF+/3qRxWLFihZw4cUI2b95sHrNnz5aLL77YzOzt16+fpKamVm6ILQggEJYCx44dk71795qH/qE4bNgwM859+/ZJWloaKRrC8qoxqFgXICAU6z8BnH80CZx//vnmdPR365ycHNH8994WvWeGFv2GHQWBcBEgwBsuV4JxIIAAAggg4EZAb/agQVx9TJw40dzRV3P2rlq1SjQo9OWXX5rHjBkzRG8QoTl7NS+Yft2MggAC4Segr1+doa9/TFpF79htBXjfeOMN+eijj+SGG26QsWPHSmIiv65bTjwjUNMCBIRq+grQPwLOCWguXf09u7S0VPTmiY888ohXja9evVq++uorU5d7ZHhFRqUQCfAbY4ig6QYBpwQ++OADmT9/vlPNlWtHvw4+YsSIcttYQQCB8BHQQE+PHj3MQ28Kob9gfvjhh/LJJ59IQUGBuWmT3ripdu3a8uCDD8rgwYPDZ/CMBIEYF9CvZU+fPl02btxYpYTm387NzTV/bG7dulWmTZvGbN4qxdiJQOgECAiFzpqeEAi2gP6+fOWVV8ry5ctlyZIlovfFuOmmm6pM16ATK5588kkzND2+d+/ewR4m7SPgtQABXq+pqIhAeAjojL1g5e7Ky8sLj5NkFAggUK1ArVq1pG/fvuZx+vRpWbBggcnLq7MQCgsL5cCBA9W28f/Zuw84qarz/+PPLr1Jr9IUC4JiiS0YiGABGxq7oAiWJL5ELIkmRiMBiSaiJkYCJtGfGoMpalQUlRALxAbYAhYsSJXepC8szP9+j7nzvzvMzE67M7O7n/N6zc7MLeee+569U5577nNY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EECgKgkwyFpVerVoKwIIIIAAAggggAACCCCAAAIIIIAAAgggEBAgwBvA4CECCCCAAAIIIIAAAggggAACCCCAAAIIIFCVBAjwVqVXi7YigAACCCCAAAIIIIAAAggggAACCCCAAAIBAQK8AQweIoAAAggggAACCCCAAAIIIIAAAggggAACVUmAAG9VerVoKwIIIIAAAggggAACCCCAAAIIIIAAAgggEBAgwBvA4CECCCCAAAIIIIAAAggggAACCCCAAAIIIFCVBAjwVqVXi7YigAACCCCAAAIIIIAAAggggAACCCCAAAIBAQK8AQweIoAAAggggAACCCCAAAIIIIAAAggggAACVUmAAG9VerVoKwIIIIAAAggggAACCCCAAAIIIIAAAgggEBCoHXjMQwQQQAABBBBAAAEEshIYP368jR07Nqs65s6da61bt86qjjBWPuqoo2zJkiXWtGlT+/TTT8PYBHUWQGDhwoXWtWvXAmyZTSKAAAIIIIAAArkRIMCbG0dqQQABBBBAAAEEEPAEtmzZYitXrszKYvfu3VmtH9bKa9ascftWVlYW1iaoN48C+l+9/fbb7c9//rMtW7Ysj1tmUwgggAACCCCAQG4FCPDm1pPaEEAAAQQQQAABBP4n0KZNG+vcuXPaHrVr8xU1bTRWSFugb9++9t5771nLli3TXpcVEEAAAQQQQACBYhLg23MxvRq0BQEEEEAAAQQQqEYCgwcPtt/85jfVZo8mTJhgW7dutTp16lSbfarJO7Ju3Tq3+yUlJTWZgX1HAAEEEEAAgWogQIC3GryI7AICCCCAAAIIIIBA+AKnnHJK+BthCwgggAACCCCAAAIIpClQmubyLI4AAggggAACCCCAAAIIIIAAAggggAACCCBQJAL04C2SF4JmIIAAAggggAACCMQXeOONN+yxxx5zMzUoVuvWre3f//63TZs2zd566y2bN2+ede/e3Y444gi7+OKL7cgjj6xQ0YsvvmjPPvusm3bmmWdaZT1x3377bXvkkUfc8qeddpqdccYZ7vEtt9xia9eutYYNG9q9997rpvl/fvzjH9vmzZvtO9/5jmvDww8/7Abv+vjjj61nz572wx/+0M477zyLTQegbWmQLy2n/WjatKn16tXL3S6//HLr0KGDv4kK99ma+JUF6/nlL3/p8tHKS76vvfaaKY3Bd7/7XTv++OPt7LPPtr322sutqoHwnnzySXvllVfs1VdftW3btrnX4Nxzz7Urr7xyj/30t+fff/DBB2593f/3v/+1evXq2WGHHeZuV1xxhbVr185fdI/7m2++2davX2/HHHOMDR8+3DT43d/+9jd788033f+Dcjjrf6F379521VVXWd26dSvUccMNN7hUG3otVfS66fVROeCAA0zzKQgggAACCCCAQJUSiFAQQAABBBBAAAEEEMiRwK9+9auI92XY3a677rqc1Prggw9G6/z8888jY8eOjT73t+Xfe/lxI17e3wrbnTNnTnT54447rsK8eE+8QGx0+VmzZkUX6dq1q5verFmz6DT/QatWrdw8LzgZueeee6Lr++3q0qVLxAuK+otHVqxYEbnkkksiXsB3j2X9dZo3bx7xAtvRdYIPsjXx6wrW8+WXX0auvfbahO0ZOHBgZNeuXZGdO3dGvPzKCZc78cQT3XL+NoL3Wt8LJEf0Ovn7GXvvBfAjzzzzTHC1Co87duzo1vWC+ZEFCxZE9t9//4R1HXXUUW6ZYAV6/WK36T/3AtnBRXmMAAIIIIAAAghUCQF68Hrf5igIIIAAAggggAACVUPgpptusqefftoNdDZgwAA7+OCDbcmSJa4n6fLly80LPtr1119vBx10kGm+yiGHHGLf+ta37N1333W9PL2goO2zzz5xd/jrr7+25557zs1Tz1svQBh3uUQTP/3002hvYy1TWlpq6u3qBXOjvVq9IKcNGjTIvOCxq0a9dvVcbdT21RN16tSprpeq1lMP13HjxiXapGViEq+yq6++2tR7Vz1gjz32WNcbWbZPPPGE7dixw1566SW79dZbTfv4z3/+07xAq3nBXGvTpo3NmDHDZs6cad4vINf796GHHnI9eWO3o97QL7zwgpvcoEEDu+iii1yv3e3btzsPL7Brq1evtrPOOsv1ktZrmah4wX7zAva2bNky84K81rdvX9cD+fXXX3d1lZeX2+zZs9021NPbL+r1qx7HkyZNsk2bNln9+vVt2LBhbrZ68FIQQAABBBBAAIEqJ1AlwtA0EgEEEEAAAQQQQKBKCAR78A4ZMiTy/vvvp3XzLpffYz+DvUy9L9sRL6gb+eijjyos5wUgI/3794/2zDz99NMrzP/9738fnacewInKn/70p+hyXlC1wmKp9OBV+3TzLvN3PUe99AGuF+7ChQujdQV7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Oj/60Y+iwV0FaV999VU3KnSilRQAfu6556xv3742d+5c27Bhg/3617+2P/7xj3usEmbde2yMCaEKnHLKKaHWT+X5EQjzdQyz7vzosJUwBfj/CFOXuhFAAAEECiVQmGQxhdpbtosAAggggAACGQmoZ60/CIhSI6RS/PQM/fr1s3bt2iVdRYOePfDAA24Z9e7QNjp16pR0Hc1s1qyZ6+XrL/jYY4+5Hsb+c92HWXdwOzxGAAEEEEAAAQQQQAABBAohQA/eQqizTQQQqHECGgzq6aeftnnz5rnbl19+6S5x33///U23Cy+80PVwrHEwRbzDSi+gyzhnzZpls2fPtr322st69+7tbuecc47FG9RDozg/+uijbq/OP/9869+/f9w91OBQCkQqADp//nzTIBlHH320q1sDmKnuG264wV2Grp6vl19+edx6/Inq6arBz9TOBQsWmEYiVlsHDBjg6vWXy+a+adOmpl5P+j+eMWOGG+Sqbdu2Catct26d/etf/3LzL7roItfTNuHC3oxJkyZFZ5999tl2+OGHR59X9kA9eOW3efNm+853vmNr166tkEIizLoraxvzkwtk8r97yy23uNdYI5ffe++9CTeg993f/va3Lm3H4sWLXQ5n5WjW/5fSqbzyyiv2j3/8w60/evRoS/b/nMn7QcKGxcy4+eabTZev6rj9/ve/7wYP/Otf/2pvvPGGzZkzx50c6dWrl40cOdKOPfbYmLW/efraa6+59xM9+9WvfmVff/21u9cgh0qvohzW6sXeuXPnCusrh/X//d//ufeODz74wOW3PtAbTFN5r0888USXD7vCCgmezJw50x588EF79913TdZ6D9tvv/2sT58+dt1117n3z9hVK3sddbWAUt7oyoEvvvjCli5d6trfo0cP0+3SSy9NeBKosrr9tih1i97n1X7l/dYAlgcddJDJWyem9F6fqPiv2zHHHGPDhw83Ddyl9/Q333zT9FlQu3bt6HuxrnioW7duoqqYnqFApsdlKv8fuf6czrStldFkeuzreNfxpfv//ve/Vq9ePXfc69i/4oorKj0p67crk2PfX1fH+F/+8hfXBv/7sY5JvXcccMABrj36TqS2UQovkMn7pf85qzEQ/KL3SV11paLvunqt45VcHDOZ/p9nelzF2w+mIZBPgRIvH14knxtkWwgggEDky1kW+eoji2xabaUd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" />

<!-- rnb-plot-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dzYXZlKFxuICBcIi4uL0ludGVybWVkaWFyeUZpbGVzL3NpbXBsaWZpZWRfcGVyZm9ybWFuY2VfbWV0cmljcnVsZV90eXBlcy5wbmdcIixcbiAgcGxvdCA9IGZpZyxcbiAgc2NhbGUgPSAxLFxuICB3aWR0aCA9IDYuNSxcbiAgaGVpZ2h0ID0gNCxcbiAgdW5pdHMgPSBjKFwiaW5cIiksXG4gIGRwaSA9IDMwMCxcbiAgbGltaXRzaXplID0gVFJVRVxuKVxuYGBgIn0= -->

```r
ggsave(
  "../IntermediaryFiles/simplified_performance_metricrule_types.png",
  plot = fig,
  scale = 1,
  width = 6.5,
  height = 4,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE
)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


visualizing change when add the tnv rule

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzX29uZV90cyA8LSBhY2N1cmFjeV9zY29yZXNbYWNjdXJhY3lfc2NvcmVzJHRlc3Rpbmdfc2V0X2luZGV4PT0xLF1cblxudG52IDwtIHJlcChOQSwgNjMpXG5cbmZvciAoaSBpbiAxOm5yb3coYWNjdXJhY3lfc2NvcmVzX29uZV90cykpIHtcbiAgdG52W2ldIDwtIGFzLm51bWVyaWMoc3Vic3RyKGFjY3VyYWN5X3Njb3Jlc19vbmVfdHMkdG9vbGNvbWJvW2ldLCA1LCA1KSlcbn1cblxubnVtX25vX3RudiA8LSA2My0oc3VtKHRudikpXG5gYGAifQ== -->

```r
accuracy_scores_one_ts <- accuracy_scores[accuracy_scores$testing_set_index==1,]

tnv <- rep(NA, 63)

for (i in 1:nrow(accuracy_scores_one_ts)) {
  tnv[i] <- as.numeric(substr(accuracy_scores_one_ts$toolcombo[i], 5, 5))
}

num_no_tnv <- 63-(sum(tnv))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubm9fdG52X3Byb3AgPC0gZGF0YS5mcmFtZSh0b29sX2NvbWJvPWFjY3VyYWN5X3Njb3Jlc19vbmVfdHMkdG9vbGNvbWJvW2dyZXAoXCIwXCIsIHRudildLFxuICAgICAgICAgICAgICAgICAgICAgICAgICBwcm9wX3ZpcmFsPXJlcCgwLCBudW1fbm9fdG52KSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgTUNDPXJlcCgwLCBudW1fbm9fdG52KSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgcHJlY2lzaW9uPXJlcCgwLCBudW1fbm9fdG52KSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgcmVjYWxsPXJlcCgwLCBudW1fbm9fdG52KSlcblxudGNfdG52IDwtIFwiMFwiXG50Y19ub190bnYgPC0gXCIwXCJcbnByb3BfdG52IDwtIDBcbmogPC0gMVxuXG5mb3IgKGkgaW4gZ3JlcChcIjBcIiwgdG52KSkge1xuICB0Y19ub190bnYgPC0gYXMuY2hhcmFjdGVyKGFjY3VyYWN5X3Njb3Jlc19vbmVfdHMkdG9vbGNvbWJvW2ldKVxuICB0Y190bnYgPC0gdGNfbm9fdG52XG4gIHN1YnN0cih0Y190bnYsIDUsIDUpIDwtIFwiMVwiXG4gIFxuICBwcm9wX3RudiA8LSBhY2N1cmFjeV9zY29yZXNfb25lX3RzJHByb3BfdmlyYWxbZ3JlcCh0Y190bnYsYWNjdXJhY3lfc2NvcmVzX29uZV90cyR0b29sY29tYm8pXVsxXVxuICBub190bnZfcHJvcCRwcm9wX3ZpcmFsW2pdIDwtXG4gICAgKHByb3BfdG52LWFjY3VyYWN5X3Njb3Jlc19vbmVfdHMkcHJvcF92aXJhbFtpXSkvcHJvcF90bnYqMTAwXG4gIFxuICBwcm9wX3RudiA8LSBhY2N1cmFjeV9zY29yZXNfb25lX3RzJE1DQ1tncmVwKHRjX3RudixhY2N1cmFjeV9zY29yZXNfb25lX3RzJHRvb2xjb21ibyldWzFdXG4gIG5vX3Rudl9wcm9wJE1DQ1tqXSA8LVxuICAgIChwcm9wX3Rudi1hY2N1cmFjeV9zY29yZXNfb25lX3RzJE1DQ1tpXSkvcHJvcF90bnYqMTAwXG4gIFxuICBwcm9wX3RudiA8LSBhY2N1cmFjeV9zY29yZXNfb25lX3RzJHByZWNpc2lvbltncmVwKHRjX3RudixhY2N1cmFjeV9zY29yZXNfb25lX3RzJHRvb2xjb21ibyldWzFdXG4gIG5vX3Rudl9wcm9wJHByZWNpc2lvbltqXSA8LVxuICAgIChwcm9wX3Rudi1hY2N1cmFjeV9zY29yZXNfb25lX3RzJHByZWNpc2lvbltpXSkvcHJvcF90bnYqMTAwXG4gIFxuICBwcm9wX3RudiA8LSBhY2N1cmFjeV9zY29yZXNfb25lX3RzJHJlY2FsbFtncmVwKHRjX3RudixhY2N1cmFjeV9zY29yZXNfb25lX3RzJHRvb2xjb21ibyldWzFdXG4gIG5vX3Rudl9wcm9wJHJlY2FsbFtqXSA8LVxuICAgIChwcm9wX3Rudi1hY2N1cmFjeV9zY29yZXNfb25lX3RzJHJlY2FsbFtpXSkvcHJvcF90bnYqMTAwXG4gIFxuICBcbiAgXG4gIGogPC0gaisxXG59XG5cblxuXG5tZWFuKG5vX3Rudl9wcm9wJHByb3BfdmlyYWwpXG5zZChub190bnZfcHJvcCRwcm9wX3ZpcmFsKVxuXG5ub190bnZfcHJvcF9tZWx0IDwtIG5vX3Rudl9wcm9wICU+JSBcbiAgc2VsZWN0KHByZWNpc2lvbiwgcmVjYWxsLCBNQ0MsIHByb3BfdmlyYWwsIHRvb2xfY29tYm8pICU+JVxuICBwaXZvdF9sb25nZXIoY29scz1jKHByZWNpc2lvbiwgcmVjYWxsLCBNQ0MsIHByb3BfdmlyYWwpLCBcbiAgICAgICAgICAgICAgIG5hbWVzX3RvPVwicGVyZm9ybWFuY2VfbWV0cmljXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJwZXJmb3JtYW5jZV9tZXRyaWNfcGVyY2VudF9kaWZmXCIpXG5cbmBgYCJ9 -->

```r
no_tnv_prop <- data.frame(tool_combo=accuracy_scores_one_ts$toolcombo[grep("0", tnv)],
                          prop_viral=rep(0, num_no_tnv),
                          MCC=rep(0, num_no_tnv),
                          precision=rep(0, num_no_tnv),
                          recall=rep(0, num_no_tnv))

tc_tnv <- "0"
tc_no_tnv <- "0"
prop_tnv <- 0
j <- 1

for (i in grep("0", tnv)) {
  tc_no_tnv <- as.character(accuracy_scores_one_ts$toolcombo[i])
  tc_tnv <- tc_no_tnv
  substr(tc_tnv, 5, 5) <- "1"
  
  prop_tnv <- accuracy_scores_one_ts$prop_viral[grep(tc_tnv,accuracy_scores_one_ts$toolcombo)][1]
  no_tnv_prop$prop_viral[j] <-
    (prop_tnv-accuracy_scores_one_ts$prop_viral[i])/prop_tnv*100
  
  prop_tnv <- accuracy_scores_one_ts$MCC[grep(tc_tnv,accuracy_scores_one_ts$toolcombo)][1]
  no_tnv_prop$MCC[j] <-
    (prop_tnv-accuracy_scores_one_ts$MCC[i])/prop_tnv*100
  
  prop_tnv <- accuracy_scores_one_ts$precision[grep(tc_tnv,accuracy_scores_one_ts$toolcombo)][1]
  no_tnv_prop$precision[j] <-
    (prop_tnv-accuracy_scores_one_ts$precision[i])/prop_tnv*100
  
  prop_tnv <- accuracy_scores_one_ts$recall[grep(tc_tnv,accuracy_scores_one_ts$toolcombo)][1]
  no_tnv_prop$recall[j] <-
    (prop_tnv-accuracy_scores_one_ts$recall[i])/prop_tnv*100
  
  
  
  j <- j+1
}



mean(no_tnv_prop$prop_viral)
sd(no_tnv_prop$prop_viral)

no_tnv_prop_melt <- no_tnv_prop %>% 
  select(precision, recall, MCC, prop_viral, tool_combo) %>%
  pivot_longer(cols=c(precision, recall, MCC, prop_viral), 
               names_to="performance_metric",
               values_to="performance_metric_percent_diff")

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


visualizing when add the tuning addition rule

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzX29uZV90cyA8LSBhY2N1cmFjeV9zY29yZXNbYWNjdXJhY3lfc2NvcmVzJHRlc3Rpbmdfc2V0X2luZGV4PT0xLF1cblxudHYgPC0gcmVwKE5BLCA2MylcblxuZm9yIChpIGluIDE6bnJvdyhhY2N1cmFjeV9zY29yZXNfb25lX3RzKSkge1xuICB0dltpXSA8LSBhcy5udW1lcmljKHN1YnN0cihhY2N1cmFjeV9zY29yZXNfb25lX3RzJHRvb2xjb21ib1tpXSwgMSwgMSkpXG59XG5cbm51bV9ub190diA8LSA2My0oc3VtKHR2KSlcblxubm9fdHZfcHJvcCA8LSBkYXRhLmZyYW1lKHRvb2xfY29tYm89YWNjdXJhY3lfc2NvcmVzX29uZV90cyR0b29sY29tYm9bZ3JlcChcIjBcIiwgdHYpXSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgcHJvcF92aXJhbD1yZXAoMCwgbnVtX25vX3R2KSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgTUNDPXJlcCgwLCBudW1fbm9fdHYpLFxuICAgICAgICAgICAgICAgICAgICAgICAgICBwcmVjaXNpb249cmVwKDAsIG51bV9ub190diksXG4gICAgICAgICAgICAgICAgICAgICAgICAgIHJlY2FsbD1yZXAoMCwgbnVtX25vX3R2KSlcblxudGNfdHYgPC0gXCIwXCJcbnRjX25vX3R2IDwtIFwiMFwiXG5wcm9wX3R2IDwtIDBcbmogPC0gMVxuXG5mb3IgKGkgaW4gZ3JlcChcIjBcIiwgdHYpKSB7XG4gIHRjX25vX3R2IDwtIGFzLmNoYXJhY3RlcihhY2N1cmFjeV9zY29yZXNfb25lX3RzJHRvb2xjb21ib1tpXSlcbiAgdGNfdHYgPC0gdGNfbm9fdHZcbiAgc3Vic3RyKHRjX3R2LCAxLCAxKSA8LSBcIjFcIlxuICBcbiAgcHJvcF90diA8LSBhY2N1cmFjeV9zY29yZXNfb25lX3RzJHByb3BfdmlyYWxbZ3JlcCh0Y190dixhY2N1cmFjeV9zY29yZXNfb25lX3RzJHRvb2xjb21ibyldWzFdXG4gIG5vX3R2X3Byb3AkcHJvcF92aXJhbFtqXSA8LVxuICAgIChwcm9wX3R2LWFjY3VyYWN5X3Njb3Jlc19vbmVfdHMkcHJvcF92aXJhbFtpXSkvcHJvcF90dioxMDBcbiAgXG4gIHByb3BfdHYgPC0gYWNjdXJhY3lfc2NvcmVzX29uZV90cyRNQ0NbZ3JlcCh0Y190dixhY2N1cmFjeV9zY29yZXNfb25lX3RzJHRvb2xjb21ibyldWzFdXG4gIG5vX3R2X3Byb3AkTUNDW2pdIDwtXG4gICAgKHByb3BfdHYtYWNjdXJhY3lfc2NvcmVzX29uZV90cyRNQ0NbaV0pL3Byb3BfdHYqMTAwXG4gIFxuICBwcm9wX3R2IDwtIGFjY3VyYWN5X3Njb3Jlc19vbmVfdHMkcHJlY2lzaW9uW2dyZXAodGNfdHYsYWNjdXJhY3lfc2NvcmVzX29uZV90cyR0b29sY29tYm8pXVsxXVxuICBub190dl9wcm9wJHByZWNpc2lvbltqXSA8LVxuICAgIChwcm9wX3R2LWFjY3VyYWN5X3Njb3Jlc19vbmVfdHMkcHJlY2lzaW9uW2ldKS9wcm9wX3R2KjEwMFxuICBcbiAgcHJvcF90diA8LSBhY2N1cmFjeV9zY29yZXNfb25lX3RzJHJlY2FsbFtncmVwKHRjX3R2LGFjY3VyYWN5X3Njb3Jlc19vbmVfdHMkdG9vbGNvbWJvKV1bMV1cbiAgbm9fdHZfcHJvcCRyZWNhbGxbal0gPC1cbiAgICBhYnMocHJvcF90di1hY2N1cmFjeV9zY29yZXNfb25lX3RzJHJlY2FsbFtpXSkvcHJvcF90dioxMDBcbiAgXG4gIFxuICBcbiAgaiA8LSBqKzFcbn1cblxubWVhbihub190dl9wcm9wJHByb3BfdmlyYWwpXG5zZChub190dl9wcm9wJHByb3BfdmlyYWwpXG5cbm5vX3R2X3Byb3BfbWVsdCA8LSBub190dl9wcm9wICU+JSBcbiAgc2VsZWN0KHByZWNpc2lvbiwgcmVjYWxsLCBNQ0MsIHByb3BfdmlyYWwsIHRvb2xfY29tYm8pICU+JVxuICBwaXZvdF9sb25nZXIoY29scz1jKHByZWNpc2lvbiwgcmVjYWxsLCBNQ0MsIHByb3BfdmlyYWwpLCBcbiAgICAgICAgICAgICAgIG5hbWVzX3RvPVwicGVyZm9ybWFuY2VfbWV0cmljXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJwZXJmb3JtYW5jZV9tZXRyaWNfcGVyY2VudF9kaWZmXCIpXG5gYGAifQ== -->

```r
accuracy_scores_one_ts <- accuracy_scores[accuracy_scores$testing_set_index==1,]

tv <- rep(NA, 63)

for (i in 1:nrow(accuracy_scores_one_ts)) {
  tv[i] <- as.numeric(substr(accuracy_scores_one_ts$toolcombo[i], 1, 1))
}

num_no_tv <- 63-(sum(tv))

no_tv_prop <- data.frame(tool_combo=accuracy_scores_one_ts$toolcombo[grep("0", tv)],
                          prop_viral=rep(0, num_no_tv),
                          MCC=rep(0, num_no_tv),
                          precision=rep(0, num_no_tv),
                          recall=rep(0, num_no_tv))

tc_tv <- "0"
tc_no_tv <- "0"
prop_tv <- 0
j <- 1

for (i in grep("0", tv)) {
  tc_no_tv <- as.character(accuracy_scores_one_ts$toolcombo[i])
  tc_tv <- tc_no_tv
  substr(tc_tv, 1, 1) <- "1"
  
  prop_tv <- accuracy_scores_one_ts$prop_viral[grep(tc_tv,accuracy_scores_one_ts$toolcombo)][1]
  no_tv_prop$prop_viral[j] <-
    (prop_tv-accuracy_scores_one_ts$prop_viral[i])/prop_tv*100
  
  prop_tv <- accuracy_scores_one_ts$MCC[grep(tc_tv,accuracy_scores_one_ts$toolcombo)][1]
  no_tv_prop$MCC[j] <-
    (prop_tv-accuracy_scores_one_ts$MCC[i])/prop_tv*100
  
  prop_tv <- accuracy_scores_one_ts$precision[grep(tc_tv,accuracy_scores_one_ts$toolcombo)][1]
  no_tv_prop$precision[j] <-
    (prop_tv-accuracy_scores_one_ts$precision[i])/prop_tv*100
  
  prop_tv <- accuracy_scores_one_ts$recall[grep(tc_tv,accuracy_scores_one_ts$toolcombo)][1]
  no_tv_prop$recall[j] <-
    abs(prop_tv-accuracy_scores_one_ts$recall[i])/prop_tv*100
  
  
  
  j <- j+1
}

mean(no_tv_prop$prop_viral)
sd(no_tv_prop$prop_viral)

no_tv_prop_melt <- no_tv_prop %>% 
  select(precision, recall, MCC, prop_viral, tool_combo) %>%
  pivot_longer(cols=c(precision, recall, MCC, prop_viral), 
               names_to="performance_metric",
               values_to="performance_metric_percent_diff")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubm9fdHZfYW5kX250dl9wcm9wX21lbHQgPC0gcmJpbmQobm9fdHZfcHJvcF9tZWx0LCBub190bnZfcHJvcF9tZWx0KVxuXG5ub190dl9hbmRfbnR2X3Byb3BfbWVsdCRzZXQgPC0gYyhyZXAoXCJ0dW5pbmcgYWRkaXRpb25cIiwgbnJvdyhub190dl9wcm9wX21lbHQpKSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJlcChcInR1bmluZyByZW1vdmFsXCIsIG5yb3cobm9fdG52X3Byb3BfbWVsdCkpKVxuXG5wYWwgPC0gZ2d0aGVtZXM6OnRhYmxlYXVfY29sb3JfcGFsKHBhbGV0dGU9XCJUYWJsZWF1IDEwXCIsIHR5cGU9XCJyZWd1bGFyXCIpXG5cblxuZmlnIDwtIGdncGxvdChub190dl9hbmRfbnR2X3Byb3BfbWVsdCwgXG4gICAgICAgICAgICAgIGFlcyh4PXBlcmZvcm1hbmNlX21ldHJpYywgeT1wZXJmb3JtYW5jZV9tZXRyaWNfcGVyY2VudF9kaWZmLFxuICAgICAgICAgICAgICAgICAgY29sb3I9c2V0LCBmaWxsPXNldCkpICtcbiAgZ2VvbV9ib3hwbG90KCkgK1xuICAjZ2VvbV9wb2ludCgpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgYXhpcy50aWNrcy55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCIsXG4gICAgYXhpcy50ZXh0Lnk9ZWxlbWVudF90ZXh0KHNpemU9MTApLFxuICAgIGF4aXMudGV4dC54PWVsZW1lbnRfdGV4dChzaXplPTEyKSxcbiAgICBsZWdlbmQudGV4dD1lbGVtZW50X3RleHQoc2l6ZT0xMiksXG4gICAgYXhpcy50aXRsZT1lbGVtZW50X3RleHQoc2l6ZT0xMiksXG4gICAgc3RyaXAuYmFja2dyb3VuZCA9IGVsZW1lbnRfcmVjdChmaWxsPVwid2hpdGVcIiwgY29sb3I9XCJncmV5XCIpLFxuICAgIHN0cmlwLnRleHQgPSBlbGVtZW50X3RleHQoY29sb3I9XCJibGFja1wiLCBzaXplPTEyKVxuICApICtcbiAgeGxhYihcIlBlcmZvcm1hbmNlIE1ldHJpY1wiKSArXG4gIHlsYWIoXCJQZXJjZW50IERpZmZlcmVuY2VcIikgK1xuICBnZW9tX2hsaW5lKHlpbnRlcmNlcHQgPSAxLCBjb2xvcj1cImdyZXlcIikgK1xuICBzY2FsZV9maWxsX21hbnVhbChuYW1lPVwiXCIsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShjKHBhbCg2KVtjKDYsMSldKSwgMC40KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhwYWwoNilbYyg2LDEpXSksIDAuNikpXG5cbmZpZ1xuYGBgIn0= -->

```r
no_tv_and_ntv_prop_melt <- rbind(no_tv_prop_melt, no_tnv_prop_melt)

no_tv_and_ntv_prop_melt$set <- c(rep("tuning addition", nrow(no_tv_prop_melt)),
                                 rep("tuning removal", nrow(no_tnv_prop_melt)))

pal <- ggthemes::tableau_color_pal(palette="Tableau 10", type="regular")


fig <- ggplot(no_tv_and_ntv_prop_melt, 
              aes(x=performance_metric, y=performance_metric_percent_diff,
                  color=set, fill=set)) +
  geom_boxplot() +
  #geom_point() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=10),
    axis.text.x=element_text(size=12),
    legend.text=element_text(size=12),
    axis.title=element_text(size=12),
    strip.background = element_rect(fill="white", color="grey"),
    strip.text = element_text(color="black", size=12)
  ) +
  xlab("Performance Metric") +
  ylab("Percent Difference") +
  geom_hline(yintercept = 1, color="grey") +
  scale_fill_manual(name="",
                     values = alpha(c(pal(6)[c(6,1)]), 0.4)) +
  scale_color_manual(name="",
                     values = alpha(c(pal(6)[c(6,1)]), 0.6))

fig
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


### Visualize performance across all the rules 


## Visualize how the precision, recall, and MCC changes across pipelines.
precision

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzIDwtIGFjY3VyYWN5X3Njb3Jlc1tvcmRlcihhY2N1cmFjeV9zY29yZXMkcHJlY2lzaW9uLCBkZWNyZWFzaW5nPUYpLF1cbmFjY3VyYWN5X3Njb3JlcyR0b29sY29tYm8gPC0gZmFjdG9yKGFjY3VyYWN5X3Njb3JlcyR0b29sY29tYm8sIGxldmVscyA9IHVuaXF1ZShhY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvKSlcblxuZ2dwbG90KGFjY3VyYWN5X3Njb3JlcywgYWVzKHg9dG9vbGNvbWJvLCB5PXByZWNpc2lvbiwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9bnVtcnVsZXMsIGZpbGw9bnVtcnVsZXMpKSArXG4gIGdlb21fcG9pbnQoYWxwaGE9MC41KSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCwgYW5nbGUgPSA5MCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTYpLFxuICApICtcbiAgeGxhYihcIlRvb2wgQ29tYmluYXRpb24gKHR2LCBEVkYsIHRudiwgVkIsIFZTLCBWUzIpXCIpICtcbiAgeWxhYihcIlByZWNpc2lvblwiKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJOdW1iZXIgb2YgUnVsZXNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMC41KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIk51bWJlciBvZiBSdWxlc1wiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAxKSlcblxuZ2dwbG90KGFjY3VyYWN5X3Njb3JlcywgYWVzKHg9dG9vbGNvbWJvLCB5PXByZWNpc2lvbiwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9cnVsZXR5cGUsIGZpbGw9cnVsZXR5cGUpKSArXG4gIGdlb21fYm94cGxvdChhbHBoYT0wLjUpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgYXhpcy50aWNrcy55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCIsXG4gICAgYXhpcy50ZXh0Lnk9ZWxlbWVudF90ZXh0KHNpemU9MTQpLFxuICAgIGF4aXMudGV4dC54PWVsZW1lbnRfdGV4dChzaXplPTE0LCBhbmdsZSA9IDkwKSxcbiAgICBsZWdlbmQudGV4dD1lbGVtZW50X3RleHQoc2l6ZT0xMiksXG4gICAgYXhpcy50aXRsZT1lbGVtZW50X3RleHQoc2l6ZT0xNiksXG4gICkgK1xuICB4bGFiKFwiVG9vbCBDb21iaW5hdGlvbiAodHYsIERWRiwgdG52LCBWQiwgVlMsIFZTMilcIikgK1xuICB5bGFiKFwiUHJlY2lzaW9uXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYobWFnbWEoNylbMzo2XSksIFwiZ3JleVwiKSwgMC40KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYobWFnbWEoNylbMzo2XSksIFwiZ3JleVwiKSwgMC42KSlcbmBgYCJ9 -->

```r
accuracy_scores <- accuracy_scores[order(accuracy_scores$precision, decreasing=F),]
accuracy_scores$toolcombo <- factor(accuracy_scores$toolcombo, levels = unique(accuracy_scores$toolcombo))

ggplot(accuracy_scores, aes(x=toolcombo, y=precision, 
                                  color=numrules, fill=numrules)) +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Tool Combination (tv, DVF, tnv, VB, VS, VS2)") +
  ylab("Precision") +
  scale_fill_manual(name="Number of Rules",
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name="Number of Rules",
                     values = alpha(rev(viridis(6)), 1))

ggplot(accuracy_scores, aes(x=toolcombo, y=precision, 
                                  color=ruletype, fill=ruletype)) +
  geom_boxplot(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Tool Combination (tv, DVF, tnv, VB, VS, VS2)") +
  ylab("Precision") +
  scale_fill_manual(name="",
                     values = alpha(c(rev(magma(7)[3:6]), "grey"), 0.4)) +
  scale_color_manual(name="",
                     values = alpha(c(rev(magma(7)[3:6]), "grey"), 0.6))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubWVhbl9wcmVjaXNpb24gPC0gYWNjdXJhY3lfc2NvcmVzICU+JSBcbiAgc2VsZWN0KHByZWNpc2lvbiwgdG9vbGNvbWJvKSAlPiVcbiAgZ3JvdXBfYnkodG9vbGNvbWJvKSAlPiVcbiAgc3VtbWFyaXNlKG1lYW4gPSBtZWFuKHByZWNpc2lvbikpXG5gYGAifQ== -->

```r
mean_precision <- accuracy_scores %>% 
  select(precision, toolcombo) %>%
  group_by(toolcombo) %>%
  summarise(mean = mean(precision))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


getting set of tools with roughly equivalent precision

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRwcmVjaXNpb25bYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIwIDAgMSAwIDAgMVwiXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJHByZWNpc2lvblthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjAgMCAxIDEgMCAwXCJdKVxuXG50LnRlc3QoYWNjdXJhY3lfc2NvcmVzJHByZWNpc2lvblthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjAgMCAxIDAgMCAxXCJdLCBcbiAgICAgICBhY2N1cmFjeV9zY29yZXMkcHJlY2lzaW9uW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMCAwIDAgMSAwIDBcIl0pXG5cbnQudGVzdChhY2N1cmFjeV9zY29yZXMkcHJlY2lzaW9uW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMCAwIDEgMCAwIDFcIl0sIFxuICAgICAgIGFjY3VyYWN5X3Njb3JlcyRwcmVjaXNpb25bYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIwIDAgMCAwIDAgMVwiXSlcblxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRwcmVjaXNpb25bYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIwIDAgMSAwIDAgMVwiXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJHByZWNpc2lvblthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjAgMCAxIDEgMCAxXCJdKVxuYGBgIn0= -->

```r
t.test(accuracy_scores$precision[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$precision[accuracy_scores$toolcombo=="0 0 1 1 0 0"])

t.test(accuracy_scores$precision[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$precision[accuracy_scores$toolcombo=="0 0 0 1 0 0"])

t.test(accuracy_scores$precision[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$precision[accuracy_scores$toolcombo=="0 0 0 0 0 1"])

t.test(accuracy_scores$precision[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$precision[accuracy_scores$toolcombo=="0 0 1 1 0 1"])
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

accounting for multiple t-tests with Bonferroni, then the following have the same precision:
- vs2+tnv, vb+tnv, vb, and vs2


recall

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzIDwtIGFjY3VyYWN5X3Njb3Jlc1tvcmRlcihhY2N1cmFjeV9zY29yZXMkcmVjYWxsLCBkZWNyZWFzaW5nPUYpLF1cbmFjY3VyYWN5X3Njb3JlcyR0b29sY29tYm8gPC0gZmFjdG9yKGFjY3VyYWN5X3Njb3JlcyR0b29sY29tYm8sIGxldmVscyA9IHVuaXF1ZShhY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvKSlcblxuZ2dwbG90KGFjY3VyYWN5X3Njb3JlcywgYWVzKHg9dG9vbGNvbWJvLCB5PXJlY2FsbCwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9bnVtcnVsZXMsIGZpbGw9bnVtcnVsZXMpKSArXG4gIGdlb21fcG9pbnQoYWxwaGE9MC41KSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCwgYW5nbGUgPSA5MCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTYpLFxuICApICtcbiAgeGxhYihcIlRvb2wgQ29tYmluYXRpb24gKHR2LCBEVkYsIHRudiwgVkIsIFZTLCBWUzIpXCIpICtcbiAgeWxhYihcIlJlY2FsbFwiKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJOdW1iZXIgb2YgUnVsZXNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMC41KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIk51bWJlciBvZiBSdWxlc1wiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAxKSlcblxuZ2dwbG90KGFjY3VyYWN5X3Njb3JlcywgYWVzKHg9dG9vbGNvbWJvLCB5PXJlY2FsbCwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9cnVsZXR5cGUsIGZpbGw9cnVsZXR5cGUpKSArXG4gIGdlb21fYm94cGxvdChhbHBoYT0wLjUpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgYXhpcy50aWNrcy55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCIsXG4gICAgYXhpcy50ZXh0Lnk9ZWxlbWVudF90ZXh0KHNpemU9MTQpLFxuICAgIGF4aXMudGV4dC54PWVsZW1lbnRfdGV4dChzaXplPTE0LCBhbmdsZSA9IDkwKSxcbiAgICBsZWdlbmQudGV4dD1lbGVtZW50X3RleHQoc2l6ZT0xMiksXG4gICAgYXhpcy50aXRsZT1lbGVtZW50X3RleHQoc2l6ZT0xNiksXG4gICkgK1xuICB4bGFiKFwiVG9vbCBDb21iaW5hdGlvbiAodHYsIERWRiwgdG52LCBWQiwgVlMsIFZTMilcIikgK1xuICB5bGFiKFwiUmVjYWxsXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYobWFnbWEoNylbMzo2XSksIFwiZ3JleVwiKSwgMC40KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYobWFnbWEoNylbMzo2XSksIFwiZ3JleVwiKSwgMC42KSlcbmBgYCJ9 -->

```r
accuracy_scores <- accuracy_scores[order(accuracy_scores$recall, decreasing=F),]
accuracy_scores$toolcombo <- factor(accuracy_scores$toolcombo, levels = unique(accuracy_scores$toolcombo))

ggplot(accuracy_scores, aes(x=toolcombo, y=recall, 
                                  color=numrules, fill=numrules)) +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Tool Combination (tv, DVF, tnv, VB, VS, VS2)") +
  ylab("Recall") +
  scale_fill_manual(name="Number of Rules",
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name="Number of Rules",
                     values = alpha(rev(viridis(6)), 1))

ggplot(accuracy_scores, aes(x=toolcombo, y=recall, 
                                  color=ruletype, fill=ruletype)) +
  geom_boxplot(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Tool Combination (tv, DVF, tnv, VB, VS, VS2)") +
  ylab("Recall") +
  scale_fill_manual(name="",
                     values = alpha(c(rev(magma(7)[3:6]), "grey"), 0.4)) +
  scale_color_manual(name="",
                     values = alpha(c(rev(magma(7)[3:6]), "grey"), 0.6))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubWVhbl9yZWNhbGwgPC0gYWNjdXJhY3lfc2NvcmVzICU+JSBcbiAgc2VsZWN0KHJlY2FsbCwgdG9vbGNvbWJvKSAlPiVcbiAgZ3JvdXBfYnkodG9vbGNvbWJvKSAlPiVcbiAgc3VtbWFyaXNlKG1lYW4gPSBtZWFuKHJlY2FsbCkpXG5gYGAifQ== -->

```r
mean_recall <- accuracy_scores %>% 
  select(recall, toolcombo) %>%
  group_by(toolcombo) %>%
  summarise(mean = mean(recall))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


getting set of tools with roughly equivalent recall

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRyZWNhbGxbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIxIDEgMCAwIDEgMVwiXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJHJlY2FsbFthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjEgMSAwIDEgMSAxXCJdKVxuXG50LnRlc3QoYWNjdXJhY3lfc2NvcmVzJHJlY2FsbFthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjEgMSAwIDAgMSAxXCJdLCBcbiAgICAgICBhY2N1cmFjeV9zY29yZXMkcmVjYWxsW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMSAxIDAgMCAwIDFcIl0pXG5cbnQudGVzdChhY2N1cmFjeV9zY29yZXMkcmVjYWxsW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMSAxIDAgMCAxIDFcIl0sIFxuICAgICAgIGFjY3VyYWN5X3Njb3JlcyRyZWNhbGxbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIxIDEgMCAxIDAgMVwiXSlcblxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRyZWNhbGxbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIxIDEgMCAwIDEgMVwiXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJHJlY2FsbFthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjEgMCAwIDAgMSAxXCJdKVxuXG50LnRlc3QoYWNjdXJhY3lfc2NvcmVzJHJlY2FsbFthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjEgMSAwIDAgMSAxXCJdLCBcbiAgICAgICBhY2N1cmFjeV9zY29yZXMkcmVjYWxsW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMSAwIDAgMSAxIDFcIl0pXG5cbnQudGVzdChhY2N1cmFjeV9zY29yZXMkcmVjYWxsW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMSAxIDAgMCAxIDFcIl0sIFxuICAgICAgIGFjY3VyYWN5X3Njb3JlcyRyZWNhbGxbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIxIDEgMSAwIDEgMVwiXSlcbmBgYCJ9 -->

```r
t.test(accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 0 1 1"], 
       accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 1 1 1"])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 0 1 1"], 
       accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 0 0 1"])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 0 1 1"], 
       accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 1 0 1"])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 0 1 1"], 
       accuracy_scores$recall[accuracy_scores$toolcombo=="1 0 0 0 1 1"])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 0 1 1"], 
       accuracy_scores$recall[accuracy_scores$toolcombo=="1 0 0 1 1 1"])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 0 0 1 1"], 
       accuracy_scores$recall[accuracy_scores$toolcombo=="1 1 1 0 1 1"])
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

accounting for multiple t-tests with Bonferroni, then the following have the same recall:
- vs2+vs+dvf+tv, vs2+vs+dvf+tv+vb, vs2+dvf+tv, vs2+vs+dvf+tv, vs2+dvf+tv+vb, 


precision vs recall

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KGFjY3VyYWN5X3Njb3JlcywgYWVzKHg9cHJlY2lzaW9uLCB5PXJlY2FsbCwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9bnVtcnVsZXMsIGZpbGw9bnVtcnVsZXMpKSArXG4gIGdlb21fcG9pbnQoYWxwaGE9MC41KSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCwgYW5nbGUgPSA5MCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MjApLFxuICApICtcbiAgeGxhYihcIlByZWNpc2lvblwiKSArXG4gIHlsYWIoXCJSZWNhbGxcIikgK1xuICBzY2FsZV9maWxsX21hbnVhbChuYW1lPVwiTnVtYmVyIG9mIFJ1bGVzXCIsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShyZXYodmlyaWRpcyg2KSksIDAuNSkpICtcbiAgc2NhbGVfY29sb3JfbWFudWFsKG5hbWU9XCJOdW1iZXIgb2YgUnVsZXNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMSkpXG5cbmdncGxvdChhY2N1cmFjeV9zY29yZXMsIGFlcyh4PXByZWNpc2lvbiwgeT1yZWNhbGwsIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGNvbG9yPW51bXRvb2xzLCBmaWxsPW51bXRvb2xzKSkgK1xuICBnZW9tX3BvaW50KGFscGhhPTAuNSkgK1xuICB0aGVtZV9saWdodCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBheGlzLnRpY2tzLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIixcbiAgICBheGlzLnRleHQueT1lbGVtZW50X3RleHQoc2l6ZT0xNCksXG4gICAgYXhpcy50ZXh0Lng9ZWxlbWVudF90ZXh0KHNpemU9MTQsIGFuZ2xlID0gOTApLFxuICAgIGxlZ2VuZC50ZXh0PWVsZW1lbnRfdGV4dChzaXplPTEyKSxcbiAgICBheGlzLnRpdGxlPWVsZW1lbnRfdGV4dChzaXplPTIwKSxcbiAgKSArXG4gIHhsYWIoXCJQcmVjaXNpb25cIikgK1xuICB5bGFiKFwiUmVjYWxsXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIk51bWJlciBvZiBUb29sc1wiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHBhbCg2KSksIDAuNSkpICtcbiAgc2NhbGVfY29sb3JfbWFudWFsKG5hbWU9XCJOdW1iZXIgb2YgVG9vbHNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldihwYWwoNikpLCAxKSlcblxuYGBgIn0= -->

```r
ggplot(accuracy_scores, aes(x=precision, y=recall, 
                                  color=numrules, fill=numrules)) +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=20),
  ) +
  xlab("Precision") +
  ylab("Recall") +
  scale_fill_manual(name="Number of Rules",
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name="Number of Rules",
                     values = alpha(rev(viridis(6)), 1))

ggplot(accuracy_scores, aes(x=precision, y=recall, 
                                  color=numtools, fill=numtools)) +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=20),
  ) +
  xlab("Precision") +
  ylab("Recall") +
  scale_fill_manual(name="Number of Tools",
                     values = alpha(rev(pal(6)), 0.5)) +
  scale_color_manual(name="Number of Tools",
                     values = alpha(rev(pal(6)), 1))

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

MCC

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubWVhbl9NQ0MgPC0gYWNjdXJhY3lfc2NvcmVzICU+JSBcbiAgc2VsZWN0KE1DQywgdG9vbGNvbWJvKSAlPiVcbiAgZ3JvdXBfYnkodG9vbGNvbWJvKSAlPiVcbiAgc3VtbWFyaXNlKG1lYW4gPSBtZWFuKE1DQykpICU+JVxuICBhcnJhbmdlKG1lYW4pXG5gYGAifQ== -->

```r
mean_MCC <- accuracy_scores %>% 
  select(MCC, toolcombo) %>%
  group_by(toolcombo) %>%
  summarise(mean = mean(MCC)) %>%
  arrange(mean)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibyA8LSBmYWN0b3IoYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibywgbGV2ZWxzID0gbWVhbl9NQ0MkdG9vbGNvbWJvKVxuXG5maWcgPC0gZ2dwbG90KGFjY3VyYWN5X3Njb3JlcywgYWVzKHg9dG9vbGNvbWJvLCB5PU1DQywgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9cnVsZXR5cGUsIGZpbGw9cnVsZXR5cGUpKSArXG4gIGdlb21fYm94cGxvdChhbHBoYT0wLjUpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgYXhpcy50aWNrcy55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCIsXG4gICAgYXhpcy50ZXh0Lnk9ZWxlbWVudF90ZXh0KHNpemU9MTQpLFxuICAgIGF4aXMudGV4dC54PWVsZW1lbnRfdGV4dChzaXplPTE0LCBhbmdsZSA9IDkwKSxcbiAgICBsZWdlbmQudGV4dD1lbGVtZW50X3RleHQoc2l6ZT0xMiksXG4gICAgYXhpcy50aXRsZT1lbGVtZW50X3RleHQoc2l6ZT0xNiksXG4gICkgK1xuICB4bGFiKFwiVG9vbCBDb21iaW5hdGlvbiAodHYsIERWRiwgdG52LCBWQiwgVlMsIFZTMilcIikgK1xuICB5bGFiKFwiTUNDXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIk51bWJlciBvZiBSdWxlc1wiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYobWFnbWEoNylbMzo2XSksIFwiZ3JleVwiKSwgMC40KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIk51bWJlciBvZiBSdWxlc1wiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYobWFnbWEoNylbMzo2XSksIFwiZ3JleVwiKSwgMC42KSlcblxuZmlnXG5cbmdnc2F2ZShcbiAgXCIuLi9JbnRlcm1lZGlhcnlGaWxlcy9hbGxfTUNDX3NjYXR0ZXJwbG90LnBuZ1wiLFxuICBwbG90ID0gZmlnLFxuICBzY2FsZSA9IDEsXG4gIHdpZHRoID0gOCxcbiAgaGVpZ2h0ID0gNCxcbiAgdW5pdHMgPSBjKFwiaW5cIiwgXCJjbVwiLCBcIm1tXCIsIFwicHhcIiksXG4gIGRwaSA9IDMwMCxcbiAgbGltaXRzaXplID0gVFJVRVxuKVxuXG5nZ3Bsb3QoYWNjdXJhY3lfc2NvcmVzLCBhZXMoeD10b29sY29tYm8sIHk9TUNDLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb2xvcj1udW10b29scywgZmlsbD1udW10b29scykpICtcbiAgZ2VvbV9wb2ludChhbHBoYT0wLjUpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgYXhpcy50aWNrcy55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCIsXG4gICAgYXhpcy50ZXh0Lnk9ZWxlbWVudF90ZXh0KHNpemU9MTQpLFxuICAgIGF4aXMudGV4dC54PWVsZW1lbnRfdGV4dChzaXplPTE0LCBhbmdsZSA9IDkwKSxcbiAgICBsZWdlbmQudGV4dD1lbGVtZW50X3RleHQoc2l6ZT0xMiksXG4gICAgYXhpcy50aXRsZT1lbGVtZW50X3RleHQoc2l6ZT0yMCksXG4gICkgK1xuICB4bGFiKFwiVG9vbCBDb21iaW5hdGlvbiAodHYsIERWRiwgdG52LCBWQiwgVlMsIFZTMilcIikgK1xuICB5bGFiKFwiTUNDXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIk51bWJlciBvZiBUb29sc1wiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAwLjUpKSArXG4gIHNjYWxlX2NvbG9yX21hbnVhbChuYW1lPVwiTnVtYmVyIG9mIFRvb2xzXCIsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShyZXYodmlyaWRpcyg2KSksIDEpKVxuYGBgIn0= -->

```r
accuracy_scores$toolcombo <- factor(accuracy_scores$toolcombo, levels = mean_MCC$toolcombo)

fig <- ggplot(accuracy_scores, aes(x=toolcombo, y=MCC, 
                                  color=ruletype, fill=ruletype)) +
  geom_boxplot(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Tool Combination (tv, DVF, tnv, VB, VS, VS2)") +
  ylab("MCC") +
  scale_fill_manual(name="Number of Rules",
                     values = alpha(c(rev(magma(7)[3:6]), "grey"), 0.4)) +
  scale_color_manual(name="Number of Rules",
                     values = alpha(c(rev(magma(7)[3:6]), "grey"), 0.6))

fig

ggsave(
  "../IntermediaryFiles/all_MCC_scatterplot.png",
  plot = fig,
  scale = 1,
  width = 8,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE
)

ggplot(accuracy_scores, aes(x=toolcombo, y=MCC, 
                                  color=numtools, fill=numtools)) +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=20),
  ) +
  xlab("Tool Combination (tv, DVF, tnv, VB, VS, VS2)") +
  ylab("MCC") +
  scale_fill_manual(name="Number of Tools",
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name="Number of Tools",
                     values = alpha(rev(viridis(6)), 1))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRyZWNhbGxbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XFwxIDEgMCAwIDEgMVxcXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJHJlY2FsbFthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cXDEgMSAwIDEgMSAxXFxdKVxuXG50LnRlc3QoYWNjdXJhY3lfc2NvcmVzJHJlY2FsbFthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cXDEgMSAwIDAgMSAxXFxdLCBcbiAgICAgICBhY2N1cmFjeV9zY29yZXMkcmVjYWxsW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVxcMSAxIDAgMCAwIDFcXF0pXG5cbnQudGVzdChhY2N1cmFjeV9zY29yZXMkcmVjYWxsW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVxcMSAxIDAgMCAxIDFcXF0sIFxuICAgICAgIGFjY3VyYWN5X3Njb3JlcyRyZWNhbGxbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XFwxIDEgMCAxIDAgMVxcXSlcblxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRyZWNhbGxbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XFwxIDEgMCAwIDEgMVxcXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJHJlY2FsbFthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cXDEgMCAwIDAgMSAxXFxdKVxuXG50LnRlc3QoYWNjdXJhY3lfc2NvcmVzJHJlY2FsbFthY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cXDEgMSAwIDAgMSAxXFxdLCBcbiAgICAgICBhY2N1cmFjeV9zY29yZXMkcmVjYWxsW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVxcMSAwIDAgMSAxIDFcXF0pXG5cbnQudGVzdChhY2N1cmFjeV9zY29yZXMkcmVjYWxsW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVxcMSAxIDAgMCAxIDFcXF0sIFxuICAgICAgIGFjY3VyYWN5X3Njb3JlcyRyZWNhbGxbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XFwxIDEgMSAwIDEgMVxcXSlcbmBgYFxuYGBgIn0= -->

```r
```r
t.test(accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 0 1 1\], 
       accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 1 1 1\])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 0 1 1\], 
       accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 0 0 1\])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 0 1 1\], 
       accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 1 0 1\])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 0 1 1\], 
       accuracy_scores$recall[accuracy_scores$toolcombo==\1 0 0 0 1 1\])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 0 1 1\], 
       accuracy_scores$recall[accuracy_scores$toolcombo==\1 0 0 1 1 1\])

t.test(accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 0 0 1 1\], 
       accuracy_scores$recall[accuracy_scores$toolcombo==\1 1 1 0 1 1\])
```
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

=78

getting set of tools with roughly equivalent MCC

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRNQ0NbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIwIDAgMSAwIDAgMVwiXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJE1DQ1thY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjAgMSAxIDAgMCAxXCJdKVxuXG50LnRlc3QoYWNjdXJhY3lfc2NvcmVzJE1DQ1thY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjAgMCAxIDAgMCAxXCJdLCBcbiAgICAgICBhY2N1cmFjeV9zY29yZXMkTUNDW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMSAwIDEgMCAwIDFcIl0pXG5cbnQudGVzdChhY2N1cmFjeV9zY29yZXMkTUNDW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMCAwIDEgMCAwIDFcIl0sIFxuICAgICAgIGFjY3VyYWN5X3Njb3JlcyRNQ0NbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIwIDAgMCAwIDAgMVwiXSlcblxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRNQ0NbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIwIDAgMSAwIDAgMVwiXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJE1DQ1thY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjEgMCAxIDEgMCAxXCJdKVxuXG50LnRlc3QoYWNjdXJhY3lfc2NvcmVzJE1DQ1thY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjAgMCAxIDAgMCAxXCJdLCBcbiAgICAgICBhY2N1cmFjeV9zY29yZXMkTUNDW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMSAwIDEgMCAwIDFcIl0pXG5cbnQudGVzdChhY2N1cmFjeV9zY29yZXMkTUNDW2FjY3VyYWN5X3Njb3JlcyR0b29sY29tYm89PVwiMCAwIDEgMCAwIDFcIl0sIFxuICAgICAgIGFjY3VyYWN5X3Njb3JlcyRNQ0NbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIxIDAgMSAwIDEgMVwiXSlcblxudC50ZXN0KGFjY3VyYWN5X3Njb3JlcyRNQ0NbYWNjdXJhY3lfc2NvcmVzJHRvb2xjb21ibz09XCIwIDAgMSAwIDAgMVwiXSwgXG4gICAgICAgYWNjdXJhY3lfc2NvcmVzJE1DQ1thY2N1cmFjeV9zY29yZXMkdG9vbGNvbWJvPT1cIjEgMSAxIDAgMCAxXCJdKVxuYGBgIn0= -->

```r
t.test(accuracy_scores$MCC[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$MCC[accuracy_scores$toolcombo=="0 1 1 0 0 1"])

t.test(accuracy_scores$MCC[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$MCC[accuracy_scores$toolcombo=="1 0 1 0 0 1"])

t.test(accuracy_scores$MCC[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$MCC[accuracy_scores$toolcombo=="0 0 0 0 0 1"])

t.test(accuracy_scores$MCC[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$MCC[accuracy_scores$toolcombo=="1 0 1 1 0 1"])

t.test(accuracy_scores$MCC[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$MCC[accuracy_scores$toolcombo=="1 0 1 0 0 1"])

t.test(accuracy_scores$MCC[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$MCC[accuracy_scores$toolcombo=="1 0 1 0 1 1"])

t.test(accuracy_scores$MCC[accuracy_scores$toolcombo=="0 0 1 0 0 1"], 
       accuracy_scores$MCC[accuracy_scores$toolcombo=="1 1 1 0 0 1"])
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

accounting for multiple t-tests with Bonferroni, then the following have the same MCC:
- vs2+tnv, vs2+tnv+dvf, vs2, vs2+tnv+tv+vb, vs2+tnv+tv+vs   


making heatmap for the x-axis label

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYW5uX2NvbCA9IGRhdGEuZnJhbWUoXG4gICAgdG52ID0gY29tYm9zX2xpc3QkdHVuZV9ub3RfdmlyYWwsXG4gICAgdHYgPSBjb21ib3NfbGlzdCR0dW5lX3ZpcmFsLFxuICAgIGR2ZiA9IGNvbWJvc19saXN0JERWRixcbiAgICB2YiA9IGNvbWJvc19saXN0JFZJQlJBTlQsXG4gICAgdnMgPSBjb21ib3NfbGlzdCRWUyxcbiAgICB2czIgPSBjb21ib3NfbGlzdCRWUzJcbilcbnJvd25hbWVzKGFubl9jb2wpIDwtIGNvbWJvc19saXN0JHRvb2xjb21ib1xuXG5hbm5fY29sb3JzID0gbGlzdChcbiAgdHYgPSBjKFwiYmxhY2tcIiwgcGFsKDYpWzNdKSxcbiAgdG52ID0gYyhcImJsYWNrXCIsIHBhbCg2KVsxXSksXG4gIGR2ZiA9IGMoXCJibGFja1wiLCBwYWwoNilbMl0pLFxuICB2YiA9IGMoXCJibGFja1wiLCBwYWwoNilbNF0pLFxuICB2cyA9IGMoXCJibGFja1wiLCBwYWwoNilbNV0pLFxuICB2czIgPSBjKFwiYmxhY2tcIiwgcGFsKDYpWzZdKVxuKVxuXG5cblxuXG52aXJhbF9zY29yZXNfMSA8LSBhY2N1cmFjeV9zY29yZXMgJT4lIGZpbHRlcih0ZXN0aW5nX3NldF9pbmRleD09NSlcbiN2aXJhbF9zY29yZXNfMSA8LSBtZWFuX01DQ1xuI3JuIDwtIHZpcmFsX3Njb3Jlc18xJHRvb2xjb21ib1xuXG52aXJhbF9zY29yZXNfMSA8LSB2aXJhbF9zY29yZXNfMSAlPiUgc2VsZWN0KGMoXCJNQ0NcIikpXG4jdmlyYWxfc2NvcmVzXzEgPC0gdmlyYWxfc2NvcmVzXzEgJT4lIHNlbGVjdChjKFwibWVhblwiKSlcbnJvd25hbWVzKHZpcmFsX3Njb3Jlc18xKSA8LSByblxuXG5waGVhdG1hcDo6cGhlYXRtYXAoYXMubWF0cml4KHZpcmFsX3Njb3Jlc18xKSxcbiAgICAgICAgICAgICAgICAgIGFubm90YXRpb25fcm93ID0gYW5uX2NvbCxcbiAgICAgICAgICAgICAgICAgIGFubm90YXRpb25fY29sb3JzID0gYW5uX2NvbG9ycyxcbiAgICAgICAgICAgICAgICAgIGJvcmRlcl9jb2xvciA9IE5BLFxuICAgICAgICAgICAgICAgICAgc2hvd19yb3duYW1lcyA9IFQsXG4gICAgICAgICAgICAgICAgICBzaG93X2NvbG5hbWVzID0gRixcbiAgICAgICAgICAgICAgICAgIGFubm90YXRpb25fbGVnZW5kID0gRixcbiAgICAgICAgICAgICAgICAgIGNsdXN0ZXJfcm93cyA9IEYsXG4gICAgICAgICAgICAgICAgICBjbHVzdGVyX2NvbHMgPSBGLFxuICAgICAgICAgICAgICAgICAgY2VsbGhlaWdodD0zLFxuICAgICAgICAgICAgICAgICAgZmlsZW5hbWUgPSBcIi4uL0ludGVybWVkaWFyeUZpbGVzL2hlYXRtYXBfTUNDX2xhYmVscy5wbmdcIlxuICAgICAgICAgICAgICAgICAgKVxuYGBgIn0= -->

```r
ann_col = data.frame(
    tnv = combos_list$tune_not_viral,
    tv = combos_list$tune_viral,
    dvf = combos_list$DVF,
    vb = combos_list$VIBRANT,
    vs = combos_list$VS,
    vs2 = combos_list$VS2
)
rownames(ann_col) <- combos_list$toolcombo

ann_colors = list(
  tv = c("black", pal(6)[3]),
  tnv = c("black", pal(6)[1]),
  dvf = c("black", pal(6)[2]),
  vb = c("black", pal(6)[4]),
  vs = c("black", pal(6)[5]),
  vs2 = c("black", pal(6)[6])
)




viral_scores_1 <- accuracy_scores %>% filter(testing_set_index==5)
#viral_scores_1 <- mean_MCC
#rn <- viral_scores_1$toolcombo

viral_scores_1 <- viral_scores_1 %>% select(c("MCC"))
#viral_scores_1 <- viral_scores_1 %>% select(c("mean"))
rownames(viral_scores_1) <- rn

pheatmap::pheatmap(as.matrix(viral_scores_1),
                  annotation_row = ann_col,
                  annotation_colors = ann_colors,
                  border_color = NA,
                  show_rownames = T,
                  show_colnames = F,
                  annotation_legend = F,
                  cluster_rows = F,
                  cluster_cols = F,
                  cellheight=3,
                  filename = "../IntermediaryFiles/heatmap_MCC_labels.png"
                  )
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZmlnIDwtIGdncGxvdChhY2N1cmFjeV9zY29yZXNfbWVsdCwgYWVzKHg9bnVtcnVsZXMsIHk9cGVyZm9ybWFuY2VfbWV0cmljX3Njb3JlLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb2xvcj1udW1ydWxlcywgZmlsbD1udW1ydWxlcykpICtcbiAgZ2VvbV9ib3hwbG90KCkgK1xuICBnZW9tX3BvaW50KGFscGhhPTAuNSkgK1xuICB0aGVtZV9saWdodCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBheGlzLnRpY2tzLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIixcbiAgICBheGlzLnRleHQueT1lbGVtZW50X3RleHQoc2l6ZT0xNCksXG4gICAgYXhpcy50ZXh0Lng9ZWxlbWVudF90ZXh0KHNpemU9MTQsIGFuZ2xlID0gOTApLFxuICAgIGxlZ2VuZC50ZXh0PWVsZW1lbnRfdGV4dChzaXplPTEyKSxcbiAgICBheGlzLnRpdGxlPWVsZW1lbnRfdGV4dChzaXplPTE2KSxcbiAgKSArXG4gIHlsYWIoXCJTY29yZVwiKSArXG4gIHhsYWIoXCJOdW1iZXIgb2YgUnVsZSBTZXRzXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIk51bWJlciBvZiBSdWxlIFNldHNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMC41KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIk51bWJlciBvZiBSdWxlIFNldHNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMSkpICtcbiAgZmFjZXRfd3JhcCh+cGVyZm9ybWFuY2VfbWV0cmljKVxuXG5maWdcblxuZ2dzYXZlKFxuICBcIi4uL0ludGVybWVkaWFyeUZpbGVzL01DQ19wcmVjaXNpb25fcmVjYWxsX2JveHBsb3RzLnBuZ1wiLFxuICBwbG90ID0gZmlnLFxuICBzY2FsZSA9IDEsXG4gIHdpZHRoID0gOCxcbiAgaGVpZ2h0ID0gNCxcbiAgdW5pdHMgPSBjKFwiaW5cIiwgXCJjbVwiLCBcIm1tXCIsIFwicHhcIiksXG4gIGRwaSA9IDMwMCxcbiAgbGltaXRzaXplID0gVFJVRVxuKVxuYGBgIn0= -->

```r
fig <- ggplot(accuracy_scores_melt, aes(x=numrules, y=performance_metric_score, 
                                  color=numrules, fill=numrules)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  ylab("Score") +
  xlab("Number of Rule Sets") +
  scale_fill_manual(name="Number of Rule Sets",
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name="Number of Rule Sets",
                     values = alpha(rev(viridis(6)), 1)) +
  facet_wrap(~performance_metric)

fig

ggsave(
  "../IntermediaryFiles/MCC_precision_recall_boxplots.png",
  plot = fig,
  scale = 1,
  width = 8,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE
)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KGFjY3VyYWN5X3Njb3Jlc19tZWx0W2FjY3VyYWN5X3Njb3Jlc19tZWx0JHBlcmZvcm1hbmNlX21ldHJpYyE9XCJNQ0NcIixdLCBhZXMoeD1udW1ydWxlcywgeT1wZXJmb3JtYW5jZV9tZXRyaWNfc2NvcmUsIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGNvbG9yPW51bXJ1bGVzLCBmaWxsPW51bXJ1bGVzKSkgK1xuICBnZW9tX2JveHBsb3QoKSArXG4gIGdlb21fcG9pbnQoYWxwaGE9MC41KSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCwgYW5nbGUgPSA5MCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTYpLFxuICApICtcbiAgeWxhYihcIlNjb3JlXCIpICtcbiAgeGxhYihcIk51bWJlciBvZiBSdWxlIFNldHNcIikgK1xuICBzY2FsZV9maWxsX21hbnVhbChuYW1lPVwiTnVtYmVyIG9mIFJ1bGUgU2V0c1wiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAwLjUpKSArXG4gIHNjYWxlX2NvbG9yX21hbnVhbChuYW1lPVwiTnVtYmVyIG9mIFJ1bGUgU2V0c1wiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAxKSkgK1xuICBmYWNldF93cmFwKH5wZXJmb3JtYW5jZV9tZXRyaWMpXG5cbmdncGxvdChhY2N1cmFjeV9zY29yZXNfbWVsdFthY2N1cmFjeV9zY29yZXNfbWVsdCRwZXJmb3JtYW5jZV9tZXRyaWMhPVwiTUNDXCIsXSwgYWVzKHg9bnVtcnVsZXMsIHk9cGVyZm9ybWFuY2VfbWV0cmljX3Njb3JlLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb2xvcj1udW10b29scywgZmlsbD1udW10b29scykpICtcbiAgZ2VvbV9ib3hwbG90KCkgK1xuICBnZW9tX3BvaW50KGFscGhhPTAuNSkgK1xuICB0aGVtZV9saWdodCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBheGlzLnRpY2tzLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIixcbiAgICBheGlzLnRleHQueT1lbGVtZW50X3RleHQoc2l6ZT0xNCksXG4gICAgYXhpcy50ZXh0Lng9ZWxlbWVudF90ZXh0KHNpemU9MTQsIGFuZ2xlID0gOTApLFxuICAgIGxlZ2VuZC50ZXh0PWVsZW1lbnRfdGV4dChzaXplPTEyKSxcbiAgICBheGlzLnRpdGxlPWVsZW1lbnRfdGV4dChzaXplPTE2KSxcbiAgKSArXG4gIHlsYWIoXCJTY29yZVwiKSArXG4gIHhsYWIoXCJOdW1iZXIgb2YgUnVsZXNcIikgK1xuICBzY2FsZV9maWxsX21hbnVhbChuYW1lPVwiTnVtYmVyIG9mIFRvb2xzXCIsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShyZXYodmlyaWRpcyg2KSksIDAuNSkpICtcbiAgc2NhbGVfY29sb3JfbWFudWFsKG5hbWU9XCJOdW1iZXIgb2YgVG9vbHNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMSkpICtcbiAgZmFjZXRfd3JhcCh+cGVyZm9ybWFuY2VfbWV0cmljKVxuXG5nZ3Bsb3QoYWNjdXJhY3lfc2NvcmVzX21lbHRbYWNjdXJhY3lfc2NvcmVzX21lbHQkcGVyZm9ybWFuY2VfbWV0cmljIT1cIk1DQ1wiLF0sIGFlcyh4PW51bXRvb2xzLCB5PXBlcmZvcm1hbmNlX21ldHJpY19zY29yZSwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9bnVtcnVsZXMsIGZpbGw9bnVtcnVsZXMpKSArXG4gIGdlb21fYm94cGxvdCgpICtcbiAgZ2VvbV9wb2ludChhbHBoYT0wLjUpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgYXhpcy50aWNrcy55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCIsXG4gICAgYXhpcy50ZXh0Lnk9ZWxlbWVudF90ZXh0KHNpemU9MTQpLFxuICAgIGF4aXMudGV4dC54PWVsZW1lbnRfdGV4dChzaXplPTE0LCBhbmdsZSA9IDkwKSxcbiAgICBsZWdlbmQudGV4dD1lbGVtZW50X3RleHQoc2l6ZT0xMiksXG4gICAgYXhpcy50aXRsZT1lbGVtZW50X3RleHQoc2l6ZT0xNiksXG4gICkgK1xuICB5bGFiKFwiU2NvcmVcIikgK1xuICB4bGFiKFwiTnVtYmVyIG9mIFRvb2xzXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIk51bWJlciBvZiBSdWxlIFNldHNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMC41KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIk51bWJlciBvZiBSdWxlIFNldHNcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMSkpICtcbiAgZmFjZXRfd3JhcCh+cGVyZm9ybWFuY2VfbWV0cmljKVxuICBcbmBgYCJ9 -->

```r
ggplot(accuracy_scores_melt[accuracy_scores_melt$performance_metric!="MCC",], aes(x=numrules, y=performance_metric_score, 
                                  color=numrules, fill=numrules)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  ylab("Score") +
  xlab("Number of Rule Sets") +
  scale_fill_manual(name="Number of Rule Sets",
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name="Number of Rule Sets",
                     values = alpha(rev(viridis(6)), 1)) +
  facet_wrap(~performance_metric)

ggplot(accuracy_scores_melt[accuracy_scores_melt$performance_metric!="MCC",], aes(x=numrules, y=performance_metric_score, 
                                  color=numtools, fill=numtools)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  ylab("Score") +
  xlab("Number of Rules") +
  scale_fill_manual(name="Number of Tools",
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name="Number of Tools",
                     values = alpha(rev(viridis(6)), 1)) +
  facet_wrap(~performance_metric)

ggplot(accuracy_scores_melt[accuracy_scores_melt$performance_metric!="MCC",], aes(x=numtools, y=performance_metric_score, 
                                  color=numrules, fill=numrules)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  ylab("Score") +
  xlab("Number of Tools") +
  scale_fill_manual(name="Number of Rule Sets",
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name="Number of Rule Sets",
                     values = alpha(rev(viridis(6)), 1)) +
  facet_wrap(~performance_metric)
  
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub25lLndheSA8LSBhb3YocHJlY2lzaW9uIH4gbnVtcnVsZXMsIGRhdGEgPSBhY2N1cmFjeV9zY29yZXMpXG5cbnN1bW1hcnkob25lLndheSlcblxuVHVrZXlIU0Qob25lLndheSlcbmBgYCJ9 -->

```r
one.way <- aov(precision ~ numrules, data = accuracy_scores)

summary(one.way)

TukeyHSD(one.way)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub25lLndheSA8LSBhb3YocmVjYWxsIH4gbnVtcnVsZXMsIGRhdGEgPSBhY2N1cmFjeV9zY29yZXMpXG5cbnN1bW1hcnkob25lLndheSlcblxuVHVrZXlIU0Qob25lLndheSlcbmBgYCJ9 -->

```r
one.way <- aov(recall ~ numrules, data = accuracy_scores)

summary(one.way)

TukeyHSD(one.way)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub25lLndheSA8LSBhb3YoTUNDIH4gbnVtcnVsZXMsIGRhdGEgPSBhY2N1cmFjeV9zY29yZXMpXG5cbnN1bW1hcnkob25lLndheSlcblxuVHVrZXlIU0Qob25lLndheSlcbmBgYCJ9 -->

```r
one.way <- aov(MCC ~ numrules, data = accuracy_scores)

summary(one.way)

TukeyHSD(one.way)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub25lLndheSA8LSBhb3YoTUNDIH4gbnVtdG9vbHMsIGRhdGEgPSBhY2N1cmFjeV9zY29yZXMpXG5cbnN1bW1hcnkob25lLndheSlcblxuVHVrZXlIU0Qob25lLndheSlcbmBgYCJ9 -->

```r
one.way <- aov(MCC ~ numtools, data = accuracy_scores)

summary(one.way)

TukeyHSD(one.way)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## differences based on genome size

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX2hpZ2hfTUNDIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxuXG52aXJ1c2VzJGNvbmZ1c2lvbl9tYXRyaXhfaGlnaF9NQ0MgPC0gXCJ0cnVlIG5lZ2F0aXZlXCJcbnZpcnVzZXMkY29uZnVzaW9uX21hdHJpeF9oaWdoX01DQ1t2aXJ1c2VzJHNlcXR5cGU9PVwidmlydXNcIiAmIHZpcnVzZXMka2VlcF9zY29yZV9oaWdoX01DQzwxXSA8LSBcImZhbHNlIG5lZ2F0aXZlXCJcbnZpcnVzZXMkY29uZnVzaW9uX21hdHJpeF9oaWdoX01DQ1t2aXJ1c2VzJHNlcXR5cGU9PVwidmlydXNcIiAmIHZpcnVzZXMka2VlcF9zY29yZV9oaWdoX01DQz49MV0gPC0gXCJ0cnVlIHBvc2l0aXZlXCJcbnZpcnVzZXMkY29uZnVzaW9uX21hdHJpeF9oaWdoX01DQ1t2aXJ1c2VzJHNlcXR5cGUhPVwidmlydXNcIiAmIHZpcnVzZXMka2VlcF9zY29yZV9oaWdoX01DQz49MV0gPC0gXCJmYWxzZSBwb3NpdGl2ZVwiXG5cbnZpcnVzZXMkc2l6ZV9jbGFzcyA8LSBcIjMtNWtiXCJcbnZpcnVzZXMkc2l6ZV9jbGFzc1t2aXJ1c2VzJGNoZWNrdl9sZW5ndGg+NTAwMF0gPC0gXCI1LTEwa2JcIlxudmlydXNlcyRzaXplX2NsYXNzW3ZpcnVzZXMkY2hlY2t2X2xlbmd0aD4xMDAwMF0gPC0gXCI+MTBrYlwiXG5cbmNvbmZ1c2lvbl9ieV90YXhhIDwtIHZpcnVzZXMgJT4lIGNvdW50KGNvbmZ1c2lvbl9tYXRyaXhfaGlnaF9NQ0MsIHNlcXR5cGUsIHNpemVfY2xhc3MsIEluZGV4KVxuXG5jb2xuYW1lcyhjb25mdXNpb25fYnlfdGF4YSkgPC0gYyhcImNvbmZ1c2lvbl9tYXRyaXhcIiwgXCJzZXF0eXBlXCIsXCJzaXplXCIsIFwiaW5kZXhcIiwgXCJjb3VudFwiKVxuXG5jb25mdXNpb25fdmlyX2NhbGxlZCA8LSBjb25mdXNpb25fYnlfdGF4YSAlPiUgZmlsdGVyKGNvbmZ1c2lvbl9tYXRyaXg9PVwidHJ1ZSBwb3NpdGl2ZVwiIHwgY29uZnVzaW9uX21hdHJpeD09XCJmYWxzZSBwb3NpdGl2ZVwiKSBcblxudHlwZV9jb3VudCA8LSB2aXJ1c2VzICU+JSBjb3VudChzZXF0eXBlLCBzaXplX2NsYXNzLCBJbmRleClcblxuY29uZnVzaW9uX3Zpcl9jYWxsZWQkcGVyX3ZpcmFsIDwtIDBcblxuZm9yIChpIGluIGMoMTpucm93KGNvbmZ1c2lvbl92aXJfY2FsbGVkKSkpIHtcbiAgY29uZnVzaW9uX3Zpcl9jYWxsZWQkcGVyX3ZpcmFsW2ldIDwtIGNvbmZ1c2lvbl92aXJfY2FsbGVkJGNvdW50W2ldL3R5cGVfY291bnQkblt0eXBlX2NvdW50JHNlcXR5cGU9PWNvbmZ1c2lvbl92aXJfY2FsbGVkJHNlcXR5cGVbaV0gJiBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR5cGVfY291bnQkSW5kZXg9PWNvbmZ1c2lvbl92aXJfY2FsbGVkJGluZGV4W2ldICZcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR5cGVfY291bnQkc2l6ZV9jbGFzcz09Y29uZnVzaW9uX3Zpcl9jYWxsZWQkc2l6ZVtpXV0qMTAwXG59XG5gYGAifQ== -->

```r
viruses$keep_score_high_MCC <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$confusion_matrix_high_MCC <- "true negative"
viruses$confusion_matrix_high_MCC[viruses$seqtype=="virus" & viruses$keep_score_high_MCC<1] <- "false negative"
viruses$confusion_matrix_high_MCC[viruses$seqtype=="virus" & viruses$keep_score_high_MCC>=1] <- "true positive"
viruses$confusion_matrix_high_MCC[viruses$seqtype!="virus" & viruses$keep_score_high_MCC>=1] <- "false positive"

viruses$size_class <- "3-5kb"
viruses$size_class[viruses$checkv_length>5000] <- "5-10kb"
viruses$size_class[viruses$checkv_length>10000] <- ">10kb"

confusion_by_taxa <- viruses %>% count(confusion_matrix_high_MCC, seqtype, size_class, Index)

colnames(confusion_by_taxa) <- c("confusion_matrix", "seqtype","size", "index", "count")

confusion_vir_called <- confusion_by_taxa %>% filter(confusion_matrix=="true positive" | confusion_matrix=="false positive") 

type_count <- viruses %>% count(seqtype, size_class, Index)

confusion_vir_called$per_viral <- 0

for (i in c(1:nrow(confusion_vir_called))) {
  confusion_vir_called$per_viral[i] <- confusion_vir_called$count[i]/type_count$n[type_count$seqtype==confusion_vir_called$seqtype[i] & 
                                                                                    type_count$Index==confusion_vir_called$index[i] &
                                                                                    type_count$size_class==confusion_vir_called$size[i]]*100
}
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY29uZnVzaW9uX3Zpcl9jYWxsZWRfc3Vic2V0IDwtIGNvbmZ1c2lvbl92aXJfY2FsbGVkW2NvbmZ1c2lvbl92aXJfY2FsbGVkJHNlcXR5cGU9PVwidmlydXNcIixdXG5cbm9uZS53YXkgPC0gYW92KHBlcl92aXJhbCB+IHNpemUsIGRhdGEgPSBjb25mdXNpb25fdmlyX2NhbGxlZF9zdWJzZXQpXG5cbnN1bW1hcnkob25lLndheSlcblxuVHVrZXlIU0Qob25lLndheSlcbmBgYCJ9 -->

```r
confusion_vir_called_subset <- confusion_vir_called[confusion_vir_called$seqtype=="virus",]

one.way <- aov(per_viral ~ size, data = confusion_vir_called_subset)

summary(one.way)

TukeyHSD(one.way)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY29uZnVzaW9uX3Zpcl9jYWxsZWQgPC0gY29uZnVzaW9uX3Zpcl9jYWxsZWQgJT4lIGdyb3VwX2J5KHNlcXR5cGUsIHNpemUpICU+JVxuICBzdW1tYXJpc2UobWVhbj1tZWFuKHBlcl92aXJhbCksIFxuICAgICAgICAgICAgc2Q9c2QocGVyX3ZpcmFsKSlcblxuY29uZnVzaW9uX3Zpcl9jYWxsZWQkc2l6ZSA8LSBmYWN0b3IoY29uZnVzaW9uX3Zpcl9jYWxsZWQkc2l6ZSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGxldmVscyA9IGMoXCIzLTVrYlwiLCBcIjUtMTBrYlwiLCBcIj4xMGtiXCIpKVxuYGBgIn0= -->

```r
confusion_vir_called <- confusion_vir_called %>% group_by(seqtype, size) %>%
  summarise(mean=mean(per_viral), 
            sd=sd(per_viral))

confusion_vir_called$size <- factor(confusion_vir_called$size,
                                    levels = c("3-5kb", "5-10kb", ">10kb"))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KGNvbmZ1c2lvbl92aXJfY2FsbGVkLCBhZXMoeT1tZWFuLCB4PXNpemUsXG4gICAgICAgICAgICAgICAgICAgZmlsbD1zZXF0eXBlLFxuICAgICAgICAgICAgICAgICAgIGNvbG9yPXNlcXR5cGUpKSArXG4gIGdlb21fYmFyKHN0YXQ9XCJpZGVudGl0eVwiLCBwb3NpdGlvbj1wb3NpdGlvbl9kb2RnZSgpKSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTYpLFxuICApICtcbiAgZ2VvbV9lcnJvcmJhcihhZXMoeW1pbj1tZWFuLXNkLCB5bWF4PW1lYW4rc2QpLCB3aWR0aD0uMixcbiAgICAgICAgICAgICAgICAgcG9zaXRpb249cG9zaXRpb25fZG9kZ2UoLjkpKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldihwYWwoNikpLCAwLjUpKSArXG4gIHNjYWxlX2NvbG9yX21hbnVhbChuYW1lPVwiXCIsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShyZXYocGFsKDYpKSwgMSkpICtcbiAgeGxhYihcIkxlbmd0aFwiKSArXG4gIHlsYWIoXCJTZXF1ZW5jZXMgQ2FsbGVkIFZpcmFsICglKVwiKSBcbmBgYCJ9 -->

```r
ggplot(confusion_vir_called, aes(y=mean, x=size,
                   fill=seqtype,
                   color=seqtype)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(.9)) +
  scale_fill_manual(name="",
                     values = alpha(rev(pal(6)), 0.5)) +
  scale_color_manual(name="",
                     values = alpha(rev(pal(6)), 1)) +
  xlab("Length") +
  ylab("Sequences Called Viral (%)") 
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxuZ2dwbG90KGFjY3VyYWN5X3Njb3Jlc19tZWx0W2FjY3VyYWN5X3Njb3Jlc19tZWx0JHBlcmZvcm1hbmNlX21ldHJpYyE9XFxNQ0NcXCxdLCBhZXMoeD1udW1ydWxlcywgeT1wZXJmb3JtYW5jZV9tZXRyaWNfc2NvcmUsIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGNvbG9yPW51bXJ1bGVzLCBmaWxsPW51bXJ1bGVzKSkgK1xuICBnZW9tX2JveHBsb3QoKSArXG4gIGdlb21fcG9pbnQoYWxwaGE9MC41KSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcXGJvdHRvbVxcLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCwgYW5nbGUgPSA5MCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTYpLFxuICApICtcbiAgeWxhYihcXFNjb3JlXFwpICtcbiAgeGxhYihcXE51bWJlciBvZiBSdWxlIFNldHNcXCkgK1xuICBzY2FsZV9maWxsX21hbnVhbChuYW1lPVxcTnVtYmVyIG9mIFJ1bGUgU2V0c1xcLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAwLjUpKSArXG4gIHNjYWxlX2NvbG9yX21hbnVhbChuYW1lPVxcTnVtYmVyIG9mIFJ1bGUgU2V0c1xcLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAxKSkgK1xuICBmYWNldF93cmFwKH5wZXJmb3JtYW5jZV9tZXRyaWMpXG5cbmdncGxvdChhY2N1cmFjeV9zY29yZXNfbWVsdFthY2N1cmFjeV9zY29yZXNfbWVsdCRwZXJmb3JtYW5jZV9tZXRyaWMhPVxcTUNDXFwsXSwgYWVzKHg9bnVtcnVsZXMsIHk9cGVyZm9ybWFuY2VfbWV0cmljX3Njb3JlLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb2xvcj1udW10b29scywgZmlsbD1udW10b29scykpICtcbiAgZ2VvbV9ib3hwbG90KCkgK1xuICBnZW9tX3BvaW50KGFscGhhPTAuNSkgK1xuICB0aGVtZV9saWdodCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBheGlzLnRpY2tzLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXFxib3R0b21cXCxcbiAgICBheGlzLnRleHQueT1lbGVtZW50X3RleHQoc2l6ZT0xNCksXG4gICAgYXhpcy50ZXh0Lng9ZWxlbWVudF90ZXh0KHNpemU9MTQsIGFuZ2xlID0gOTApLFxuICAgIGxlZ2VuZC50ZXh0PWVsZW1lbnRfdGV4dChzaXplPTEyKSxcbiAgICBheGlzLnRpdGxlPWVsZW1lbnRfdGV4dChzaXplPTE2KSxcbiAgKSArXG4gIHlsYWIoXFxTY29yZVxcKSArXG4gIHhsYWIoXFxOdW1iZXIgb2YgUnVsZXNcXCkgK1xuICBzY2FsZV9maWxsX21hbnVhbChuYW1lPVxcTnVtYmVyIG9mIFRvb2xzXFwsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShyZXYodmlyaWRpcyg2KSksIDAuNSkpICtcbiAgc2NhbGVfY29sb3JfbWFudWFsKG5hbWU9XFxOdW1iZXIgb2YgVG9vbHNcXCxcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMSkpICtcbiAgZmFjZXRfd3JhcCh+cGVyZm9ybWFuY2VfbWV0cmljKVxuXG5nZ3Bsb3QoYWNjdXJhY3lfc2NvcmVzX21lbHRbYWNjdXJhY3lfc2NvcmVzX21lbHQkcGVyZm9ybWFuY2VfbWV0cmljIT1cXE1DQ1xcLF0sIGFlcyh4PW51bXRvb2xzLCB5PXBlcmZvcm1hbmNlX21ldHJpY19zY29yZSwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9bnVtcnVsZXMsIGZpbGw9bnVtcnVsZXMpKSArXG4gIGdlb21fYm94cGxvdCgpICtcbiAgZ2VvbV9wb2ludChhbHBoYT0wLjUpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgYXhpcy50aWNrcy55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFxcYm90dG9tXFwsXG4gICAgYXhpcy50ZXh0Lnk9ZWxlbWVudF90ZXh0KHNpemU9MTQpLFxuICAgIGF4aXMudGV4dC54PWVsZW1lbnRfdGV4dChzaXplPTE0LCBhbmdsZSA9IDkwKSxcbiAgICBsZWdlbmQudGV4dD1lbGVtZW50X3RleHQoc2l6ZT0xMiksXG4gICAgYXhpcy50aXRsZT1lbGVtZW50X3RleHQoc2l6ZT0xNiksXG4gICkgK1xuICB5bGFiKFxcU2NvcmVcXCkgK1xuICB4bGFiKFxcTnVtYmVyIG9mIFRvb2xzXFwpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cXE51bWJlciBvZiBSdWxlIFNldHNcXCxcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMC41KSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cXE51bWJlciBvZiBSdWxlIFNldHNcXCxcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMSkpICtcbiAgZmFjZXRfd3JhcCh+cGVyZm9ybWFuY2VfbWV0cmljKVxuICBcbmBgYFxuYGBgIn0= -->

```r
```r
ggplot(accuracy_scores_melt[accuracy_scores_melt$performance_metric!=\MCC\,], aes(x=numrules, y=performance_metric_score, 
                                  color=numrules, fill=numrules)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = \bottom\,
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  ylab(\Score\) +
  xlab(\Number of Rule Sets\) +
  scale_fill_manual(name=\Number of Rule Sets\,
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name=\Number of Rule Sets\,
                     values = alpha(rev(viridis(6)), 1)) +
  facet_wrap(~performance_metric)

ggplot(accuracy_scores_melt[accuracy_scores_melt$performance_metric!=\MCC\,], aes(x=numrules, y=performance_metric_score, 
                                  color=numtools, fill=numtools)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = \bottom\,
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  ylab(\Score\) +
  xlab(\Number of Rules\) +
  scale_fill_manual(name=\Number of Tools\,
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name=\Number of Tools\,
                     values = alpha(rev(viridis(6)), 1)) +
  facet_wrap(~performance_metric)

ggplot(accuracy_scores_melt[accuracy_scores_melt$performance_metric!=\MCC\,], aes(x=numtools, y=performance_metric_score, 
                                  color=numrules, fill=numrules)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = \bottom\,
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  ylab(\Score\) +
  xlab(\Number of Tools\) +
  scale_fill_manual(name=\Number of Rule Sets\,
                     values = alpha(rev(viridis(6)), 0.5)) +
  scale_color_manual(name=\Number of Rule Sets\,
                     values = alpha(rev(viridis(6)), 1)) +
  facet_wrap(~performance_metric)
  
```
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



## visualizing a select set of tools

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX2hpZ2hfYWxsIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBUKVxudmlydXNlcyRrZWVwX3Njb3JlX2hpZ2hfcmVjYWxsIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBUKVxudmlydXNlcyRrZWVwX3Njb3JlX2hpZ2hfcHJlY2lzaW9uIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxudmlydXNlcyRrZWVwX3Njb3JlX2hpZ2hfTUNDIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxuXG52aXJ1c2VzJHRydWVfdmlydXMgPC0gXCJub3RcIlxudmlydXNlcyR0cnVlX3ZpcnVzW3ZpcnVzZXMkc2VxdHlwZT09XCJ2aXJ1c1wiXSA8LSBcInZpcnVzXCJcblxudmlydXNlc19sb25nX3Njb3JlcyA8LSB2aXJ1c2VzICU+JSBcbiAgc2VsZWN0KGNvbnRhaW5zKFwia2VlcF9zY29yZV9oaWdoXCIpLCB0cnVlX3ZpcnVzKSAlPiVcbiAgcGl2b3RfbG9uZ2VyKGNvbHM9Y29udGFpbnMoXCJrZWVwX3Njb3JlX1wiKSwgXG4gICAgICAgICAgICAgICBuYW1lc190bz1cInJ1bGVfY29tYmluYXRpb25cIixcbiAgICAgICAgICAgICAgIHZhbHVlc190bz1cInZpcmFsX3Njb3JlXCIpICU+JSBcbiAgbXV0YXRlKHZpcmFsX3Njb3JlPWFzLmZhY3Rvcihyb3VuZCh2aXJhbF9zY29yZSkpKSAlPiVcbiAgZ3JvdXBfYnkocnVsZV9jb21iaW5hdGlvbiwgdmlyYWxfc2NvcmUsIHRydWVfdmlydXMpICU+JVxuICBzdW1tYXJpc2UobiA9IG4oKSlcblxucGFsX3B1cnBsZSA8LSBSQ29sb3JCcmV3ZXI6OmJyZXdlci5wYWwoNiwgXCJQdXJwbGVzXCIpXG5wYWxfb3JhbmdlIDwtIFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg0LCBcIk9yYW5nZXNcIilcblxuZ2dwbG90KHZpcnVzZXNfbG9uZ19zY29yZXMsIGFlcyh5PW4sIHg9cnVsZV9jb21iaW5hdGlvbixcbiAgICAgICAgICAgICAgICAgICBmaWxsPXZpcmFsX3Njb3JlKSkgK1xuICBnZW9tX2JhcihzdGF0PVwiaWRlbnRpdHlcIiwgY29sb3I9XCJibGFja1wiKSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICBjb29yZF9mbGlwKCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCJcbiAgKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKGMocmV2KHBhbF9vcmFuZ2UpLCBwYWxfcHVycGxlKSwgMC44KSkgK1xuICB4bGFiKFwiXCIpICtcbiAgeWxhYihcIk51bWJlciBvZiBTZXF1ZW5jZXNcIikgKyBcbiAgZmFjZXRfZ3JpZCh+dHJ1ZV92aXJ1cywgc2NhbGVzID0gXCJmcmVlXCIpXG5gYGAifQ== -->

```r
viruses$keep_score_high_all <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)
viruses$keep_score_high_recall <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = F,
                                              include_virsorter = T)
viruses$keep_score_high_precision <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)
viruses$keep_score_high_MCC <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$true_virus <- "not"
viruses$true_virus[viruses$seqtype=="virus"] <- "virus"

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_high"), true_virus) %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, true_virus) %>%
  summarise(n = n())

pal_purple <- RColorBrewer::brewer.pal(6, "Purples")
pal_orange <- RColorBrewer::brewer.pal(4, "Oranges")

ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_manual(name="",
                     values = alpha(c(rev(pal_orange), pal_purple), 0.8)) +
  xlab("") +
  ylab("Number of Sequences") + 
  facet_grid(~true_virus, scales = "free")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlc19sb25nX3Njb3JlcyA8LSB2aXJ1c2VzICU+JSBcbiAgc2VsZWN0KGNvbnRhaW5zKFwia2VlcF9zY29yZV9oaWdoXCIpLCBtYXhfc2NvcmVfZ3JvdXAsIHRydWVfdmlydXMpICU+JVxuICBmaWx0ZXIodHJ1ZV92aXJ1cz09XCJ2aXJ1c1wiKSAlPiVcbiAgcGl2b3RfbG9uZ2VyKGNvbHM9Y29udGFpbnMoXCJrZWVwX3Njb3JlX1wiKSwgXG4gICAgICAgICAgICAgICBuYW1lc190bz1cInJ1bGVfY29tYmluYXRpb25cIixcbiAgICAgICAgICAgICAgIHZhbHVlc190bz1cInZpcmFsX3Njb3JlXCIpICU+JSBcbiAgbXV0YXRlKHZpcmFsX3Njb3JlPWFzLmZhY3Rvcihyb3VuZCh2aXJhbF9zY29yZSkpKSAlPiVcbiAgZ3JvdXBfYnkocnVsZV9jb21iaW5hdGlvbiwgdmlyYWxfc2NvcmUsIG1heF9zY29yZV9ncm91cCkgJT4lXG4gIHN1bW1hcmlzZShuID0gbigpKVxuXG5wYWxfcHVycGxlIDwtIFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg2LCBcIlB1cnBsZXNcIilcbnBhbF9vcmFuZ2UgPC0gUkNvbG9yQnJld2VyOjpicmV3ZXIucGFsKDMsIFwiT3Jhbmdlc1wiKVxuXG5nZ3Bsb3QodmlydXNlc19sb25nX3Njb3JlcywgYWVzKHk9biwgeD1ydWxlX2NvbWJpbmF0aW9uLFxuICAgICAgICAgICAgICAgICAgIGZpbGw9dmlyYWxfc2NvcmUpKSArXG4gIGdlb21fYmFyKHN0YXQ9XCJpZGVudGl0eVwiLCBjb2xvcj1cImJsYWNrXCIpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIGNvb3JkX2ZsaXAoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIlxuICApICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYocGFsX29yYW5nZSksIHBhbF9wdXJwbGUpLCAwLjgpKSArXG4gIHhsYWIoXCJcIikgK1xuICB5bGFiKFwiTnVtYmVyIG9mIFNlcXVlbmNlc1wiKSArIFxuICBmYWNldF9ncmlkKH5tYXhfc2NvcmVfZ3JvdXAsIHNjYWxlcyA9IFwiZnJlZVwiKSArXG4gIHNjYWxlX3lfY29udGludW91cyhicmVha3MgPSBzY2FsZXM6OnByZXR0eV9icmVha3MobiA9IDIpKVxuXG52aXJ1c2VzX2xvbmdfc2NvcmVzIDwtIHZpcnVzZXMgJT4lIFxuICBzZWxlY3QoY29udGFpbnMoXCJrZWVwX3Njb3JlX2hpZ2hcIiksIG1heF9zY29yZV9ncm91cCwgdHJ1ZV92aXJ1cykgJT4lXG4gIGZpbHRlcih0cnVlX3ZpcnVzPT1cIm5vdFwiKSAlPiVcbiAgcGl2b3RfbG9uZ2VyKGNvbHM9Y29udGFpbnMoXCJrZWVwX3Njb3JlX1wiKSwgXG4gICAgICAgICAgICAgICBuYW1lc190bz1cInJ1bGVfY29tYmluYXRpb25cIixcbiAgICAgICAgICAgICAgIHZhbHVlc190bz1cInZpcmFsX3Njb3JlXCIpICU+JSBcbiAgbXV0YXRlKHZpcmFsX3Njb3JlPWFzLmZhY3Rvcihyb3VuZCh2aXJhbF9zY29yZSkpKSAlPiVcbiAgZ3JvdXBfYnkocnVsZV9jb21iaW5hdGlvbiwgdmlyYWxfc2NvcmUsIG1heF9zY29yZV9ncm91cCkgJT4lXG4gIHN1bW1hcmlzZShuID0gbigpKVxuXG5wYWxfcHVycGxlIDwtIFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg2LCBcIlB1cnBsZXNcIilcbnBhbF9vcmFuZ2UgPC0gUkNvbG9yQnJld2VyOjpicmV3ZXIucGFsKDQsIFwiT3Jhbmdlc1wiKVxuXG5nZ3Bsb3QodmlydXNlc19sb25nX3Njb3JlcywgYWVzKHk9biwgeD1ydWxlX2NvbWJpbmF0aW9uLFxuICAgICAgICAgICAgICAgICAgIGZpbGw9dmlyYWxfc2NvcmUpKSArXG4gIGdlb21fYmFyKHN0YXQ9XCJpZGVudGl0eVwiLCBjb2xvcj1cImJsYWNrXCIpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIGNvb3JkX2ZsaXAoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIlxuICApICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYocGFsX29yYW5nZSksIHBhbF9wdXJwbGUpLCAwLjgpKSArXG4gIHhsYWIoXCJcIikgK1xuICB5bGFiKFwiTnVtYmVyIG9mIFNlcXVlbmNlc1wiKSArIFxuICBmYWNldF9ncmlkKH5tYXhfc2NvcmVfZ3JvdXAsIHNjYWxlcyA9IFwiZnJlZVwiKSArXG4gIHNjYWxlX3lfY29udGludW91cyhicmVha3MgPSBzY2FsZXM6OnByZXR0eV9icmVha3MobiA9IDIpKVxuYGBgIn0= -->

```r
viruses_long_scores <- viruses %>% 
  select(contains("keep_score_high"), max_score_group, true_virus) %>%
  filter(true_virus=="virus") %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, max_score_group) %>%
  summarise(n = n())

pal_purple <- RColorBrewer::brewer.pal(6, "Purples")
pal_orange <- RColorBrewer::brewer.pal(3, "Oranges")

ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_manual(name="",
                     values = alpha(c(rev(pal_orange), pal_purple), 0.8)) +
  xlab("") +
  ylab("Number of Sequences") + 
  facet_grid(~max_score_group, scales = "free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2))

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_high"), max_score_group, true_virus) %>%
  filter(true_virus=="not") %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, max_score_group) %>%
  summarise(n = n())

pal_purple <- RColorBrewer::brewer.pal(6, "Purples")
pal_orange <- RColorBrewer::brewer.pal(4, "Oranges")

ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_manual(name="",
                     values = alpha(c(rev(pal_orange), pal_purple), 0.8)) +
  xlab("") +
  ylab("Number of Sequences") + 
  facet_grid(~max_score_group, scales = "free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





### Effect of combining tools

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX2R2ZiA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcblxudmlydXNlcyRrZWVwX3Njb3JlX2R2Zl92cyA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2R2Zl92c192YiA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2R2Zl92c192Yl92czIgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9kdmZfdnNfdmJfdnMyX3R2IDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBUKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfZHZmX3ZzX3ZiX3ZzMl90dl90bnYgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMkdHJ1ZV92aXJ1cyA8LSBcIm5vdFwiXG52aXJ1c2VzJHRydWVfdmlydXNbdmlydXNlcyRzZXF0eXBlPT1cInZpcnVzXCJdIDwtIFwidmlydXNcIlxuXG52aXJ1c2VzX2xvbmdfc2NvcmVzIDwtIHZpcnVzZXMgJT4lIFxuICBzZWxlY3QoY29udGFpbnMoXCJrZWVwX3Njb3JlX2R2ZlwiKSwgc2l6ZV9jbGFzcywgc2VxdHlwZSkgJT4lXG4gIHBpdm90X2xvbmdlcihjb2xzPWNvbnRhaW5zKFwia2VlcF9zY29yZV9cIiksIFxuICAgICAgICAgICAgICAgbmFtZXNfdG89XCJydWxlX2NvbWJpbmF0aW9uXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJ2aXJhbF9zY29yZVwiKSAlPiUgXG4gIG11dGF0ZSh2aXJhbF9zY29yZT1hcy5mYWN0b3Iocm91bmQodmlyYWxfc2NvcmUpKSkgJT4lXG4gIGdyb3VwX2J5KHJ1bGVfY29tYmluYXRpb24sIHZpcmFsX3Njb3JlLCBzaXplX2NsYXNzLCBzZXF0eXBlKSAlPiVcbiAgc3VtbWFyaXNlKG4gPSBuKCkpXG5cbnZpcnVzZXNfbG9uZ19zY29yZXMkc2l6ZV9jbGFzcyA8LSBmYWN0b3IodmlydXNlc19sb25nX3Njb3JlcyRzaXplX2NsYXNzLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbGV2ZWxzID0gYyhcIjMtNWtiXCIsIFwiNS0xMGtiXCIsIFwiPjEwa2JcIikpXG5gYGAifQ== -->

```r
viruses$keep_score_dvf <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = F,
                                              include_virsorter2 = F,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_dvf_vs <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = F,
                                              include_virsorter2 = F,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = T)

viruses$keep_score_dvf_vs_vb <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = F,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = T)

viruses$keep_score_dvf_vs_vb_vs2 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = T)

viruses$keep_score_dvf_vs_vb_vs2_tv <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = F,
                                              include_virsorter = T)

viruses$keep_score_dvf_vs_vb_vs2_tv_tnv <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$true_virus <- "not"
viruses$true_virus[viruses$seqtype=="virus"] <- "virus"

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_dvf"), size_class, seqtype) %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, size_class, seqtype) %>%
  summarise(n = n())

viruses_long_scores$size_class <- factor(viruses_long_scores$size_class,
                                    levels = c("3-5kb", "5-10kb", ">10kb"))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KHZpcnVzZXNfbG9uZ19zY29yZXMsIGFlcyh5PW4sIHg9cnVsZV9jb21iaW5hdGlvbixcbiAgICAgICAgICAgICAgICAgICBmaWxsPXZpcmFsX3Njb3JlKSkgK1xuICBnZW9tX2JhcihzdGF0PVwiaWRlbnRpdHlcIikgK1xuICB0aGVtZV9saWdodCgpICtcbiAgY29vcmRfZmxpcCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiXG4gICkgK1xuICBzY2FsZV9maWxsX2JyZXdlcihwYWxldHRlID0gXCJQdU9yXCIsICkgK1xuICB4bGFiKFwiXCIpICtcbiAgeWxhYihcIk51bWJlciBvZiBTZXF1ZW5jZXNcIikgKyBcbiAgc2NhbGVfeF9kaXNjcmV0ZShsYWJlbHM9YyhcIkRWRlwiLCBcIkRWRitWU1wiLCBcIkRWRitWUytWQlwiLCBcIkRWRitWUytWQitWUzJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcIkRWRitWUytWQitWUzIrYWRkaXRpb25cIiwgXCJEVkYrVlMrVkIrVlMyK2FkZGl0aW9uLXJlbW92YWxcIikpICtcbiAgZmFjZXRfZ3JpZCh+c2VxdHlwZSwgc2NhbGVzID0gXCJmcmVlXCIpXG5gYGAifQ== -->

```r
ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_brewer(palette = "PuOr", ) +
  xlab("") +
  ylab("Number of Sequences") + 
  scale_x_discrete(labels=c("DVF", "DVF+VS", "DVF+VS+VB", "DVF+VS+VB+VS2",
                            "DVF+VS+VB+VS2+addition", "DVF+VS+VB+VS2+addition-removal")) +
  facet_grid(~seqtype, scales = "free")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



### Effect of tuning rules - subset

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfYWxsIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgdHZfMT1ULCB0dl8yPVQsIHR2XzM9RiwgdHZfND1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBudHZfMT1ULCBudHZfMj1ULCBudHZfMz1ULCBudHZfND1GLCBudHZfNT1GLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBUKVxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfbmVpdGhlciA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfMCA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR2XzE9VCwgdHZfMj1ULCB0dl8zPUYsIHR2XzQ9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBudHZfMT1GLCBudHZfMj1GLCBudHZfMz1GLCBudHZfND1GLCBudHZfNT1GLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfMSA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR2XzE9VCwgdHZfMj1ULCB0dl8zPUYsIHR2XzQ9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9VCwgbnR2XzI9RiwgbnR2XzM9RiwgbnR2XzQ9RiwgbnR2XzU9RixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfMiA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR2XzE9VCwgdHZfMj1ULCB0dl8zPUYsIHR2XzQ9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9RiwgbnR2XzI9VCwgbnR2XzM9RiwgbnR2XzQ9RiwgbnR2XzU9RixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfMyA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR2XzE9VCwgdHZfMj1ULCB0dl8zPUYsIHR2XzQ9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9RiwgbnR2XzI9RiwgbnR2XzM9VCwgbnR2XzQ9RiwgbnR2XzU9RixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF90dl8wIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9VCwgbnR2XzI9VCwgbnR2XzM9VCwgbnR2XzQ9RiwgbnR2XzU9RixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9hbGxfdHZfMSA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPVQsIHR2XzI9RiwgdHZfMz1GLCB0dl80PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9VCwgbnR2XzI9VCwgbnR2XzM9VCwgbnR2XzQ9RiwgbnR2XzU9RixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF90dl8yIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR2XzE9RiwgdHZfMj1ULCB0dl8zPUYsIHR2XzQ9RixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBudHZfMT1ULCBudHZfMj1ULCBudHZfMz1ULCBudHZfND1GLCBudHZfNT1GLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBUKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfYWxsX3R2XzQgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgdHZfMT1GLCB0dl8yPUYsIHR2XzM9RiwgdHZfND1ULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG50dl8xPVQsIG50dl8yPVQsIG50dl8zPVQsIG50dl80PUYsIG50dl81PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cblxudmlydXNlcyR0cnVlX3ZpcnVzIDwtIFwibm90XCJcbnZpcnVzZXMkdHJ1ZV92aXJ1c1t2aXJ1c2VzJHNlcXR5cGU9PVwidmlydXNcIl0gPC0gXCJ2aXJ1c1wiXG5cbnZpcnVzZXNfbG9uZ19zY29yZXMgPC0gdmlydXNlcyAlPiUgXG4gIHNlbGVjdChjb250YWlucyhcImtlZXBfc2NvcmVfYWxsX1wiKSwgc2VxdHlwZSkgJT4lXG4gIHBpdm90X2xvbmdlcihjb2xzPWNvbnRhaW5zKFwia2VlcF9zY29yZV9cIiksIFxuICAgICAgICAgICAgICAgbmFtZXNfdG89XCJydWxlX2NvbWJpbmF0aW9uXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJ2aXJhbF9zY29yZVwiKSAlPiUgXG4gIG11dGF0ZSh2aXJhbF9zY29yZT1hcy5mYWN0b3Iocm91bmQodmlyYWxfc2NvcmUpKSkgJT4lXG4gIGdyb3VwX2J5KHJ1bGVfY29tYmluYXRpb24sIHZpcmFsX3Njb3JlLCBzZXF0eXBlKSAlPiVcbiAgc3VtbWFyaXNlKG4gPSBuKCkpIFxuXG52aXJ1c2VzX2xvbmdfc2NvcmVzJHR1bmluZ190eXBlIDwtIFwidHVuaW5nIGFkZGl0aW9uXCJcbnZpcnVzZXNfbG9uZ19zY29yZXMkdHVuaW5nX3R5cGVbZ3JlcGwoXCJudHZcIiwgdmlydXNlc19sb25nX3Njb3JlcyRydWxlX2NvbWJpbmF0aW9uKV0gPC0gXCJ0dW5pbmcgcmVtb3ZhbFwiXG52aXJ1c2VzX2xvbmdfc2NvcmVzJHR1bmluZ190eXBlW3ZpcnVzZXNfbG9uZ19zY29yZXMkcnVsZV9jb21iaW5hdGlvbj09XCJrZWVwX3Njb3JlX2FsbF9udHZfYWxsXCJdIDwtIFwiYWxsXCJcbnZpcnVzZXNfbG9uZ19zY29yZXMkdHVuaW5nX3R5cGVbdmlydXNlc19sb25nX3Njb3JlcyRydWxlX2NvbWJpbmF0aW9uPT1cImtlZXBfc2NvcmVfYWxsX250dl9uZWl0aGVyXCJdIDwtIFwibmVpdGhlclwiXG5cbmZpZyA8LSBnZ3Bsb3QodmlydXNlc19sb25nX3Njb3JlcywgYWVzKHk9biwgeD1ydWxlX2NvbWJpbmF0aW9uLFxuICAgICAgICAgICAgICAgICAgIGZpbGw9dmlyYWxfc2NvcmUpKSArXG4gIGdlb21fYmFyKHN0YXQ9XCJpZGVudGl0eVwiLCBjb2xvcj1cImJsYWNrXCIpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIGNvb3JkX2ZsaXAoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIixcbiAgICBzdHJpcC5iYWNrZ3JvdW5kID0gZWxlbWVudF9yZWN0KGZpbGw9XCJ3aGl0ZVwiLCBjb2xvcj1cImdyZXlcIiksXG4gICAgc3RyaXAudGV4dCA9IGVsZW1lbnRfdGV4dChjb2xvcj1cImJsYWNrXCIpXG4gICkgK1xuICBzY2FsZV9maWxsX21hbnVhbChuYW1lPVwiVmlyYWwgU2NvcmVcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKGMocmV2KFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg5LFwiT3Jhbmdlc1wiKVs0OjddKSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBSQ29sb3JCcmV3ZXI6OmJyZXdlci5wYWwoOSxcIlB1cnBsZXNcIilbNDo5XSBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICApLCAwLjgpKSArXG4gIHhsYWIoXCJSdWxlc2V0XCIpICtcbiAgeWxhYihcIlByb3BvcnRpb25cIikgKyBcbiAgZmFjZXRfZ3JpZCh0dW5pbmdfdHlwZX5zZXF0eXBlLCBzY2FsZXMgPSBcImZyZWVcIikgXG5cbmZpZ1xuXG5nZ3NhdmUoXG4gIFwiLi4vSW50ZXJtZWRpYXJ5RmlsZXMvdHVuaW5nX3J1bGVzZXRzLnBuZ1wiLFxuICBwbG90ID0gZmlnLFxuICBzY2FsZSA9IDEsXG4gIHdpZHRoID0gOCxcbiAgaGVpZ2h0ID0gNixcbiAgdW5pdHMgPSBjKFwiaW5cIiksXG4gIGRwaSA9IDMwMCxcbiAgbGltaXRzaXplID0gVFJVRVxuKVxuXG5gYGAifQ== -->

```r
viruses$keep_score_all_ntv_all <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              tv_1=T, tv_2=T, tv_3=F, tv_4=T,
                                              include_tuning_viral = T,
                                              ntv_1=T, ntv_2=T, ntv_3=T, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)
viruses$keep_score_all_ntv_neither <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_0 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              tv_1=T, tv_2=T, tv_3=F, tv_4=T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = F,
                                              ntv_1=F, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_1 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              tv_1=T, tv_2=T, tv_3=F, tv_4=T,
                                              include_tuning_viral = T,
                                              ntv_1=T, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_2 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              tv_1=T, tv_2=T, tv_3=F, tv_4=T,
                                              include_tuning_viral = T,
                                              ntv_1=F, ntv_2=T, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_3 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              tv_1=T, tv_2=T, tv_3=F, tv_4=T,
                                              include_tuning_viral = T,
                                              ntv_1=F, ntv_2=F, ntv_3=T, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_tv_0 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = T,
                                              ntv_1=T, ntv_2=T, ntv_3=T, ntv_4=F, ntv_5=F,
                                              include_virsorter = T)

viruses$keep_score_all_tv_1 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=T, tv_2=F, tv_3=F, tv_4=F,
                                              ntv_1=T, ntv_2=T, ntv_3=T, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_tv_2 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=F, tv_2=T, tv_3=F, tv_4=F,
                                              ntv_1=T, ntv_2=T, ntv_3=T, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_tv_4 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=F, tv_2=F, tv_3=F, tv_4=T,
                                              ntv_1=T, ntv_2=T, ntv_3=T, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)


viruses$true_virus <- "not"
viruses$true_virus[viruses$seqtype=="virus"] <- "virus"

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_all_"), seqtype) %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, seqtype) %>%
  summarise(n = n()) 

viruses_long_scores$tuning_type <- "tuning addition"
viruses_long_scores$tuning_type[grepl("ntv", viruses_long_scores$rule_combination)] <- "tuning removal"
viruses_long_scores$tuning_type[viruses_long_scores$rule_combination=="keep_score_all_ntv_all"] <- "all"
viruses_long_scores$tuning_type[viruses_long_scores$rule_combination=="keep_score_all_ntv_neither"] <- "neither"

fig <- ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill="white", color="grey"),
    strip.text = element_text(color="black")
  ) +
  scale_fill_manual(name="Viral Score",
                     values = alpha(c(rev(RColorBrewer::brewer.pal(9,"Oranges")[4:7]),
                            RColorBrewer::brewer.pal(9,"Purples")[4:9] 
                            ), 0.8)) +
  xlab("Ruleset") +
  ylab("Proportion") + 
  facet_grid(tuning_type~seqtype, scales = "free") 

fig

ggsave(
  "../IntermediaryFiles/tuning_rulesets.png",
  plot = fig,
  scale = 1,
  width = 8,
  height = 6,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE
)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->






### Effect of tuning addition rules

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF90dl9hbGwgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9hbGxfdHZfMSA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPVQsIHR2XzI9RiwgdHZfMz1GLCB0dl80PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9hbGxfdHZfMiA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPVQsIHR2XzI9VCwgdHZfMz1GLCB0dl80PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9hbGxfdHZfMyA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPVQsIHR2XzI9VCwgdHZfMz1ULCB0dl80PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9hbGxfdHZfMCA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPUYsIHR2XzI9RiwgdHZfMz1GLCB0dl80PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMkdHJ1ZV92aXJ1cyA8LSBcIm5vdFwiXG52aXJ1c2VzJHRydWVfdmlydXNbdmlydXNlcyRzZXF0eXBlPT1cInZpcnVzXCJdIDwtIFwidmlydXNcIlxuXG52aXJ1c2VzX2xvbmdfc2NvcmVzIDwtIHZpcnVzZXMgJT4lIFxuICBzZWxlY3QoY29udGFpbnMoXCJrZWVwX3Njb3JlX2FsbF90dlwiKSwgc2VxdHlwZSkgJT4lXG4gIHBpdm90X2xvbmdlcihjb2xzPWNvbnRhaW5zKFwia2VlcF9zY29yZV9cIiksIFxuICAgICAgICAgICAgICAgbmFtZXNfdG89XCJydWxlX2NvbWJpbmF0aW9uXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJ2aXJhbF9zY29yZVwiKSAlPiUgXG4gIG11dGF0ZSh2aXJhbF9zY29yZT1hcy5mYWN0b3Iocm91bmQodmlyYWxfc2NvcmUpKSkgJT4lXG4gIGdyb3VwX2J5KHJ1bGVfY29tYmluYXRpb24sIHZpcmFsX3Njb3JlLCBzZXF0eXBlKSAlPiVcbiAgc3VtbWFyaXNlKG4gPSBuKCkpXG5cbmZpZyA8LSBnZ3Bsb3QodmlydXNlc19sb25nX3Njb3JlcywgYWVzKHk9biwgeD1ydWxlX2NvbWJpbmF0aW9uLFxuICAgICAgICAgICAgICAgICAgIGZpbGw9dmlyYWxfc2NvcmUpKSArXG4gIGdlb21fYmFyKHN0YXQ9XCJpZGVudGl0eVwiLCBjb2xvcj1cImJsYWNrXCIpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIGNvb3JkX2ZsaXAoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIixcbiAgICBzdHJpcC5iYWNrZ3JvdW5kID0gZWxlbWVudF9yZWN0KGZpbGw9XCJ3aGl0ZVwiLCBjb2xvcj1cImdyZXlcIiksXG4gICAgc3RyaXAudGV4dCA9IGVsZW1lbnRfdGV4dChjb2xvcj1cImJsYWNrXCIpXG4gICkgK1xuICBzY2FsZV9maWxsX21hbnVhbChuYW1lPVwiVmlyYWwgU2NvcmVcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKGMocmV2KFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg5LFwiT3Jhbmdlc1wiKVs0OjddKSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBSQ29sb3JCcmV3ZXI6OmJyZXdlci5wYWwoOSxcIlB1cnBsZXNcIilbNDo5XSBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICApLCAwLjgpKSArXG4gIHhsYWIoXCJSdWxlc2V0XCIpICtcbiAgeWxhYihcIlByb3BvcnRpb25cIikgKyBcbiAgZmFjZXRfZ3JpZCh+c2VxdHlwZSwgc2NhbGVzID0gXCJmcmVlXCIpIFxuXG5maWdcblxuZ2dzYXZlKFxuICBcIi4uL0ludGVybWVkaWFyeUZpbGVzL3R1bmluZ19hZGRpdGlvbl9ydWxlc2V0c19hbGwucG5nXCIsXG4gIHBsb3QgPSBmaWcsXG4gIHNjYWxlID0gMSxcbiAgd2lkdGggPSA4LFxuICBoZWlnaHQgPSA0LFxuICB1bml0cyA9IGMoXCJpblwiKSxcbiAgZHBpID0gMzAwLFxuICBsaW1pdHNpemUgPSBUUlVFXG4pXG5gYGAifQ== -->

```r
viruses$keep_score_all_tv_all <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_tv_1 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=T, tv_2=F, tv_3=F, tv_4=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_tv_2 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=T, tv_2=T, tv_3=F, tv_4=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_tv_3 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=T, tv_2=T, tv_3=T, tv_4=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_tv_0 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              tv_1=F, tv_2=F, tv_3=F, tv_4=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$true_virus <- "not"
viruses$true_virus[viruses$seqtype=="virus"] <- "virus"

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_all_tv"), seqtype) %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, seqtype) %>%
  summarise(n = n())

fig <- ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill="white", color="grey"),
    strip.text = element_text(color="black")
  ) +
  scale_fill_manual(name="Viral Score",
                     values = alpha(c(rev(RColorBrewer::brewer.pal(9,"Oranges")[4:7]),
                            RColorBrewer::brewer.pal(9,"Purples")[4:9] 
                            ), 0.8)) +
  xlab("Ruleset") +
  ylab("Proportion") + 
  facet_grid(~seqtype, scales = "free") 

fig

ggsave(
  "../IntermediaryFiles/tuning_addition_rulesets_all.png",
  plot = fig,
  scale = 1,
  width = 8,
  height = 4,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE
)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


same idea but VS2 base

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX3ZzMl90dl9hbGwgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfdHZfNF90dl8yIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR2XzE9RiwgdHZfMj1ULCB0dl8zPUYsIHR2XzQ9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcblxudmlydXNlcyRrZWVwX3Njb3JlX3ZzMl90dl80X3R2XzJfdHZfMSA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPVQsIHR2XzI9VCwgdHZfMz1GLCB0dl80PVQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfdHZfMSA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPVQsIHR2XzI9RiwgdHZfMz1GLCB0dl80PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfdHZfMiA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPUYsIHR2XzI9VCwgdHZfMz1GLCB0dl80PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfdHZfMyA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPUYsIHR2XzI9RiwgdHZfMz1ULCB0dl80PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfdHZfNCA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB0dl8xPUYsIHR2XzI9RiwgdHZfMz1GLCB0dl80PVQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfdHZfMCA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcbiB2aXJ1c2VzJGtlZXBfc2NvcmVfdnMyX3R2X250diA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRikgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcblxudmlydXNlcyR0cnVlX3ZpcnVzIDwtIFwibm90XCJcbnZpcnVzZXMkdHJ1ZV92aXJ1c1t2aXJ1c2VzJHNlcXR5cGU9PVwidmlydXNcIl0gPC0gXCJ2aXJ1c1wiXG5cbnZpcnVzZXNfbG9uZ19zY29yZXMgPC0gdmlydXNlcyAlPiUgXG4gIHNlbGVjdChjb250YWlucyhcImtlZXBfc2NvcmVfdnMyX3R2X1wiKSwgc2VxdHlwZSkgJT4lXG4gIHBpdm90X2xvbmdlcihjb2xzPWNvbnRhaW5zKFwia2VlcF9zY29yZV9cIiksIFxuICAgICAgICAgICAgICAgbmFtZXNfdG89XCJydWxlX2NvbWJpbmF0aW9uXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJ2aXJhbF9zY29yZVwiKSAlPiUgXG4gIG11dGF0ZSh2aXJhbF9zY29yZT1hcy5mYWN0b3Iocm91bmQodmlyYWxfc2NvcmUpKSkgJT4lXG4gIGdyb3VwX2J5KHJ1bGVfY29tYmluYXRpb24sIHZpcmFsX3Njb3JlLCBzZXF0eXBlKSAlPiVcbiAgc3VtbWFyaXNlKG4gPSBuKCkpXG5gYGAifQ== -->

```r
viruses$keep_score_vs2_tv_all <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_tv_4_tv_2 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=F, tv_2=T, tv_3=F, tv_4=T,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_tv_4_tv_2_tv_1 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=T, tv_2=T, tv_3=F, tv_4=T,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_tv_1 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=T, tv_2=F, tv_3=F, tv_4=F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_tv_2 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=F, tv_2=T, tv_3=F, tv_4=F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_tv_3 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=F, tv_2=F, tv_3=T, tv_4=F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_tv_4 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=F, tv_2=F, tv_3=F, tv_4=T,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_tv_0 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)
                                              
 viruses$keep_score_vs2_tv_ntv <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)                                             

viruses$true_virus <- "not"
viruses$true_virus[viruses$seqtype=="virus"] <- "virus"

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_vs2_tv_"), seqtype) %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, seqtype) %>%
  summarise(n = n())
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZmlnIDwtIGdncGxvdCh2aXJ1c2VzX2xvbmdfc2NvcmVzLCBhZXMoeT1uLCB4PXJ1bGVfY29tYmluYXRpb24sXG4gICAgICAgICAgICAgICAgICAgZmlsbD12aXJhbF9zY29yZSkpICtcbiAgZ2VvbV9iYXIoc3RhdD1cImlkZW50aXR5XCIsIGNvbG9yPVwiYmxhY2tcIikgK1xuICB0aGVtZV9saWdodCgpICtcbiAgY29vcmRfZmxpcCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIHN0cmlwLmJhY2tncm91bmQgPSBlbGVtZW50X3JlY3QoZmlsbD1cIndoaXRlXCIsIGNvbG9yPVwiZ3JleVwiKSxcbiAgICBzdHJpcC50ZXh0ID0gZWxlbWVudF90ZXh0KGNvbG9yPVwiYmxhY2tcIilcbiAgKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJWaXJhbCBTY29yZVwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYoUkNvbG9yQnJld2VyOjpicmV3ZXIucGFsKDksXCJPcmFuZ2VzXCIpWzQ6N10pLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg5LFwiUHVycGxlc1wiKVs0OjldIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICksIDAuOCkpICtcbiAgeGxhYihcIlJ1bGVzZXRcIikgK1xuICB5bGFiKFwiUHJvcG9ydGlvblwiKSArIFxuICBmYWNldF9ncmlkKH5zZXF0eXBlLCBzY2FsZXMgPSBcImZyZWVcIikgXG5cbmZpZ1xuYGBgIn0= -->

```r
fig <- ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill="white", color="grey"),
    strip.text = element_text(color="black")
  ) +
  scale_fill_manual(name="Viral Score",
                     values = alpha(c(rev(RColorBrewer::brewer.pal(9,"Oranges")[4:7]),
                            RColorBrewer::brewer.pal(9,"Purples")[4:9] 
                            ), 0.8)) +
  xlab("Ruleset") +
  ylab("Proportion") + 
  facet_grid(~seqtype, scales = "free") 

fig
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



### Effect of tuning removal rules

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfYWxsIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBUKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfYWxsX250dl8xIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG50dl8xPVQsIG50dl8yPUYsIG50dl8zPUYsIG50dl80PUYsIG50dl81PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9hbGxfbnR2XzIgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9RiwgbnR2XzI9VCwgbnR2XzM9RiwgbnR2XzQ9RiwgbnR2XzU9RixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfMyA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBudHZfMT1GLCBudHZfMj1GLCBudHZfMz1ULCBudHZfND1GLCBudHZfNT1GLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBUKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfYWxsX250dl80IDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG50dl8xPUYsIG50dl8yPUYsIG50dl8zPUYsIG50dl80PVQsIG50dl81PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IFQpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9hbGxfbnR2XzUgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9RiwgbnR2XzI9RiwgbnR2XzM9RiwgbnR2XzQ9RiwgbnR2XzU9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfMCA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBudHZfMT1GLCBudHZfMj1GLCBudHZfMz1GLCBudHZfND1GLCBudHZfNT1GLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBUKVxuXG52aXJ1c2VzJHRydWVfdmlydXMgPC0gXCJub3RcIlxudmlydXNlcyR0cnVlX3ZpcnVzW3ZpcnVzZXMkc2VxdHlwZT09XCJ2aXJ1c1wiXSA8LSBcInZpcnVzXCJcblxudmlydXNlc19sb25nX3Njb3JlcyA8LSB2aXJ1c2VzICU+JSBcbiAgc2VsZWN0KGNvbnRhaW5zKFwia2VlcF9zY29yZV9hbGxfbnR2XCIpLCBzZXF0eXBlKSAlPiVcbiAgcGl2b3RfbG9uZ2VyKGNvbHM9Y29udGFpbnMoXCJrZWVwX3Njb3JlX1wiKSwgXG4gICAgICAgICAgICAgICBuYW1lc190bz1cInJ1bGVfY29tYmluYXRpb25cIixcbiAgICAgICAgICAgICAgIHZhbHVlc190bz1cInZpcmFsX3Njb3JlXCIpICU+JSBcbiAgbXV0YXRlKHZpcmFsX3Njb3JlPWFzLmZhY3Rvcihyb3VuZCh2aXJhbF9zY29yZSkpKSAlPiVcbiAgZ3JvdXBfYnkocnVsZV9jb21iaW5hdGlvbiwgdmlyYWxfc2NvcmUsIHNlcXR5cGUpICU+JVxuICBzdW1tYXJpc2UobiA9IG4oKSlcblxuZmlnIDwtIGdncGxvdCh2aXJ1c2VzX2xvbmdfc2NvcmVzLCBhZXMoeT1uLCB4PXJ1bGVfY29tYmluYXRpb24sXG4gICAgICAgICAgICAgICAgICAgZmlsbD12aXJhbF9zY29yZSkpICtcbiAgZ2VvbV9iYXIoc3RhdD1cImlkZW50aXR5XCIsIGNvbG9yPVwiYmxhY2tcIikgK1xuICB0aGVtZV9saWdodCgpICtcbiAgY29vcmRfZmxpcCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIHN0cmlwLmJhY2tncm91bmQgPSBlbGVtZW50X3JlY3QoZmlsbD1cIndoaXRlXCIsIGNvbG9yPVwiZ3JleVwiKSxcbiAgICBzdHJpcC50ZXh0ID0gZWxlbWVudF90ZXh0KGNvbG9yPVwiYmxhY2tcIilcbiAgKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJWaXJhbCBTY29yZVwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEoYyhyZXYoUkNvbG9yQnJld2VyOjpicmV3ZXIucGFsKDksXCJPcmFuZ2VzXCIpWzQ6N10pLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg5LFwiUHVycGxlc1wiKVs0OjldIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICksIDAuOCkpICtcbiAgeGxhYihcIlJ1bGVzZXRcIikgK1xuICB5bGFiKFwiUHJvcG9ydGlvblwiKSArIFxuICBmYWNldF9ncmlkKH5zZXF0eXBlLCBzY2FsZXMgPSBcImZyZWVcIikgXG5cbmZpZ1xuXG5nZ3NhdmUoXG4gIFwiLi4vSW50ZXJtZWRpYXJ5RmlsZXMvdHVuaW5nX3JlbW92YWxfcnVsZXNldHMtYWxsLnBuZ1wiLFxuICBwbG90ID0gZmlnLFxuICBzY2FsZSA9IDEsXG4gIHdpZHRoID0gOCxcbiAgaGVpZ2h0ID0gNCxcbiAgdW5pdHMgPSBjKFwiaW5cIiksXG4gIGRwaSA9IDMwMCxcbiAgbGltaXRzaXplID0gVFJVRVxuKVxuYGBgIn0= -->

```r
viruses$keep_score_all_ntv_all <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_1 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              ntv_1=T, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_2 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              ntv_1=F, ntv_2=T, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_3 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              ntv_1=F, ntv_2=F, ntv_3=T, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_4 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              ntv_1=F, ntv_2=F, ntv_3=F, ntv_4=T, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_5 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              ntv_1=F, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = T)

viruses$keep_score_all_ntv_0 <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = T,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              ntv_1=F, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = T)

viruses$true_virus <- "not"
viruses$true_virus[viruses$seqtype=="virus"] <- "virus"

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_all_ntv"), seqtype) %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, seqtype) %>%
  summarise(n = n())

fig <- ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill="white", color="grey"),
    strip.text = element_text(color="black")
  ) +
  scale_fill_manual(name="Viral Score",
                     values = alpha(c(rev(RColorBrewer::brewer.pal(9,"Oranges")[4:7]),
                            RColorBrewer::brewer.pal(9,"Purples")[4:9] 
                            ), 0.8)) +
  xlab("Ruleset") +
  ylab("Proportion") + 
  facet_grid(~seqtype, scales = "free") 

fig

ggsave(
  "../IntermediaryFiles/tuning_removal_rulesets-all.png",
  plot = fig,
  scale = 1,
  width = 8,
  height = 4,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE
)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfMCwgdmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfMSlcbmBgYCJ9 -->

```r
table(viruses$keep_score_all_ntv_0, viruses$keep_score_all_ntv_1)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfNCwgdmlydXNlcyRrZWVwX3Njb3JlX2FsbF9udHZfNSlcbmBgYCJ9 -->

```r
table(viruses$keep_score_all_ntv_4, viruses$keep_score_all_ntv_5)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KHZpcnVzZXNfbG9uZ19zY29yZXMsIGFlcyh5PW4sIHg9cnVsZV9jb21iaW5hdGlvbixcbiAgICAgICAgICAgICAgICAgICBmaWxsPXZpcmFsX3Njb3JlKSkgK1xuICBnZW9tX2JhcihzdGF0PVwiaWRlbnRpdHlcIiwgY29sb3I9XCJibGFja1wiKSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICBjb29yZF9mbGlwKCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCJcbiAgKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKGMocmV2KHBhbF9vcmFuZ2UpLCBwYWxfcHVycGxlKSwgMSkpICtcbiAgeGxhYihcIlwiKSArXG4gIHlsYWIoXCJOdW1iZXIgb2YgU2VxdWVuY2VzXCIpICsgXG4gIyBzY2FsZV94X2Rpc2NyZXRlKGxhYmVscz1jKFwiRFZGXCIsIFwiRFZGK1ZTXCIsIFwiRFZGK1ZTK1ZCXCIsIFwiRFZGK1ZTK1ZCK1ZTMlwiLFxuICMgICAgICAgICAgICAgICAgICAgICAgICAgIFwiRFZGK1ZTK1ZCK1ZTMithZGRpdGlvblwiLCBcIkRWRitWUytWQitWUzIrYWRkaXRpb24tcmVtb3ZhbFwiKSkgKyBcbiAgZmFjZXRfZ3JpZChzaXplX2NsYXNzfnNlcXR5cGUsIHNjYWxlcyA9IFwiZnJlZVwiKVxuYGBgIn0= -->

```r
ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_manual(name="",
                     values = alpha(c(rev(pal_orange), pal_purple), 1)) +
  xlab("") +
  ylab("Number of Sequences") + 
 # scale_x_discrete(labels=c("DVF", "DVF+VS", "DVF+VS+VB", "DVF+VS+VB+VS2",
 #                          "DVF+VS+VB+VS2+addition", "DVF+VS+VB+VS2+addition-removal")) + 
  facet_grid(size_class~seqtype, scales = "free")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





same idea but with a VS2 base

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX3ZzMl9udHZfYWxsIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfdnMyX250dl8xIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG50dl8xPVQsIG50dl8yPUYsIG50dl8zPUYsIG50dl80PUYsIG50dl81PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfbnR2XzIgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9RiwgbnR2XzI9VCwgbnR2XzM9RiwgbnR2XzQ9RiwgbnR2XzU9RixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcblxudmlydXNlcyRrZWVwX3Njb3JlX3ZzMl9udHZfMyA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBudHZfMT1GLCBudHZfMj1GLCBudHZfMz1ULCBudHZfND1GLCBudHZfNT1GLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfdnMyX250dl80IDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG50dl8xPUYsIG50dl8yPUYsIG50dl8zPUYsIG50dl80PVQsIG50dl81PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfbnR2XzUgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbnR2XzE9RiwgbnR2XzI9RiwgbnR2XzM9RiwgbnR2XzQ9RiwgbnR2XzU9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcblxudmlydXNlcyRrZWVwX3Njb3JlX3ZzMl9udHZfMCA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBudHZfMT1GLCBudHZfMj1GLCBudHZfMz1GLCBudHZfND1GLCBudHZfNT1GLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfdnMyX250dl8wX3R2IDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG50dl8xPUYsIG50dl8yPUYsIG50dl8zPUYsIG50dl80PUYsIG50dl81PUYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV92czJfbnR2X2FsbF90diA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcblxudmlydXNlcyRrZWVwX3Njb3JlX3ZzMl9udHZfYWxsX3R2XzRfMl8xIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHR2XzE9VCwgdHZfMj1ULCB0dl8zPUYsIHR2XzQ9VCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcblxudmlydXNlcyR0cnVlX3ZpcnVzIDwtIFwibm90XCJcbnZpcnVzZXMkdHJ1ZV92aXJ1c1t2aXJ1c2VzJHNlcXR5cGU9PVwidmlydXNcIl0gPC0gXCJ2aXJ1c1wiXG5cbnZpcnVzZXNfbG9uZ19zY29yZXMgPC0gdmlydXNlcyAlPiUgXG4gIHNlbGVjdChjb250YWlucyhcImtlZXBfc2NvcmVfdnMyX250dlwiKSwgc2VxdHlwZSkgJT4lXG4gIHBpdm90X2xvbmdlcihjb2xzPWNvbnRhaW5zKFwia2VlcF9zY29yZV92czJcIiksIFxuICAgICAgICAgICAgICAgbmFtZXNfdG89XCJydWxlX2NvbWJpbmF0aW9uXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJ2aXJhbF9zY29yZVwiKSAlPiUgXG4gIG11dGF0ZSh2aXJhbF9zY29yZT1hcy5mYWN0b3Iocm91bmQodmlyYWxfc2NvcmUpKSkgJT4lXG4gIGdyb3VwX2J5KHJ1bGVfY29tYmluYXRpb24sIHZpcmFsX3Njb3JlLCBzZXF0eXBlKSAlPiVcbiAgc3VtbWFyaXNlKG4gPSBuKCkpXG5gYGAifQ== -->

```r
viruses$keep_score_vs2_ntv_all <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_1 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              ntv_1=T, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_2 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              ntv_1=F, ntv_2=T, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_3 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              ntv_1=F, ntv_2=F, ntv_3=T, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_4 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              ntv_1=F, ntv_2=F, ntv_3=F, ntv_4=T, ntv_5=F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_5 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              ntv_1=F, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_0 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              ntv_1=F, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_0_tv <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              ntv_1=F, ntv_2=F, ntv_3=F, ntv_4=F, ntv_5=F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_all_tv <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$keep_score_vs2_ntv_all_tv_4_2_1 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              tv_1=T, tv_2=T, tv_3=F, tv_4=T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$true_virus <- "not"
viruses$true_virus[viruses$seqtype=="virus"] <- "virus"

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_vs2_ntv"), seqtype) %>%
  pivot_longer(cols=contains("keep_score_vs2"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, seqtype) %>%
  summarise(n = n())
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KHZpcnVzZXNfbG9uZ19zY29yZXMsIGFlcyh5PW4sIHg9cnVsZV9jb21iaW5hdGlvbixcbiAgICAgICAgICAgICAgICAgICBmaWxsPXZpcmFsX3Njb3JlKSkgK1xuICBnZW9tX2JhcihzdGF0PVwiaWRlbnRpdHlcIiwgY29sb3I9XCJibGFja1wiKSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICBjb29yZF9mbGlwKCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCJcbiAgKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKGMocmV2KHBhbF9vcmFuZ2UpLCBwYWxfcHVycGxlKSwgMC43KSkgK1xuICB4bGFiKFwiXCIpICtcbiAgeWxhYihcIk51bWJlciBvZiBTZXF1ZW5jZXNcIikgKyBcbiAgZmFjZXRfZ3JpZCh+c2VxdHlwZSwgc2NhbGVzID0gXCJmcmVlXCIpXG5gYGAifQ== -->

```r
ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_manual(name="",
                     values = alpha(c(rev(pal_orange), pal_purple), 0.7)) +
  xlab("") +
  ylab("Number of Sequences") + 
  facet_grid(~seqtype, scales = "free")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




comparing each tool's score

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlcyRrZWVwX3Njb3JlX2luZF92czIgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9pbmRfdHYgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX2RlZXB2aXJmaW5kZXIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlicmFudCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ192aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IEYpXG5cbnZpcnVzZXMka2VlcF9zY29yZV9pbmRfZHZmIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfaW5kX3ZiIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlcywgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxuXG52aXJ1c2VzJGtlZXBfc2NvcmVfaW5kX250diA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcblxudmlydXNlcyRrZWVwX3Njb3JlX2luZF92cyA8LSBnZXR0aW5nX3ZpcmFsX3NldF8xKHZpcnVzZXMsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gVClcblxudmlydXNlcyR0cnVlX3ZpcnVzIDwtIFwibm90XCJcbnZpcnVzZXMkdHJ1ZV92aXJ1c1t2aXJ1c2VzJHNlcXR5cGU9PVwidmlydXNcIl0gPC0gXCJ2aXJ1c1wiXG5cbnZpcnVzZXNfbG9uZ19zY29yZXMgPC0gdmlydXNlcyAlPiUgXG4gIHNlbGVjdChjb250YWlucyhcImtlZXBfc2NvcmVfaW5kX1wiKSwgc2VxdHlwZSkgJT4lXG4gIHBpdm90X2xvbmdlcihjb2xzPWNvbnRhaW5zKFwia2VlcF9zY29yZV9cIiksIFxuICAgICAgICAgICAgICAgbmFtZXNfdG89XCJydWxlX2NvbWJpbmF0aW9uXCIsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJ2aXJhbF9zY29yZVwiKSAlPiUgXG4gIG11dGF0ZSh2aXJhbF9zY29yZT1hcy5mYWN0b3Iocm91bmQodmlyYWxfc2NvcmUpKSkgJT4lXG4gIGdyb3VwX2J5KHJ1bGVfY29tYmluYXRpb24sIHZpcmFsX3Njb3JlLCBzZXF0eXBlKSAlPiVcbiAgc3VtbWFyaXNlKG4gPSBuKCkpXG5gYGAifQ== -->

```r
viruses$keep_score_ind_vs2 <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_ind_tv <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_ind_dvf <- getting_viral_set_1(viruses, include_deepvirfinder = T,
                                              include_vibrant = F,
                                              include_virsorter2 = F,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_ind_vb <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = T,
                                              include_virsorter2 = F,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = F)

viruses$keep_score_ind_ntv <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = F,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses$keep_score_ind_vs <- getting_viral_set_1(viruses, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = F,
                                              include_tuning_viral = F,
                                              include_tuning_not_viral = F,
                                              include_virsorter = T)

viruses$true_virus <- "not"
viruses$true_virus[viruses$seqtype=="virus"] <- "virus"

viruses_long_scores <- viruses %>% 
  select(contains("keep_score_ind_"), seqtype) %>%
  pivot_longer(cols=contains("keep_score_"), 
               names_to="rule_combination",
               values_to="viral_score") %>% 
  mutate(viral_score=as.factor(round(viral_score))) %>%
  group_by(rule_combination, viral_score, seqtype) %>%
  summarise(n = n())
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

Fig 4


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodmlydXNlcyR0cnVlX3ZpcnVzLCB2aXJ1c2VzJGtlZXBfc2NvcmVfaW5kX2R2ZilcbnRhYmxlKHZpcnVzZXMkdHJ1ZV92aXJ1cywgdmlydXNlcyRrZWVwX3Njb3JlX2luZF92czIpXG50YWJsZSh2aXJ1c2VzJHRydWVfdmlydXMsIHZpcnVzZXMka2VlcF9zY29yZV9pbmRfdnMpXG50YWJsZSh2aXJ1c2VzJHRydWVfdmlydXMsIHZpcnVzZXMka2VlcF9zY29yZV9pbmRfdmIpXG5gYGAifQ== -->

```r
table(viruses$true_virus, viruses$keep_score_ind_dvf)
table(viruses$true_virus, viruses$keep_score_ind_vs2)
table(viruses$true_virus, viruses$keep_score_ind_vs)
table(viruses$true_virus, viruses$keep_score_ind_vb)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

seems like false positive rate is way high for 0.5 across all


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodmlydXNlcyR0cnVlX3ZpcnVzLCB2aXJ1c2VzJGtlZXBfc2NvcmVfaW5kX2R2Zit2aXJ1c2VzJGtlZXBfc2NvcmVfaW5kX3ZzMilcbmBgYCJ9 -->

```r
table(viruses$true_virus, viruses$keep_score_ind_dvf+viruses$keep_score_ind_vs2)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodmlydXNlcyR0cnVlX3ZpcnVzLCB2aXJ1c2VzJGtlZXBfc2NvcmVfaW5kX2R2Zit2aXJ1c2VzJGtlZXBfc2NvcmVfaW5kX3ZzMit2aXJ1c2VzJGtlZXBfc2NvcmVfaW5kX250dilcbmBgYCJ9 -->

```r
table(viruses$true_virus, viruses$keep_score_ind_dvf+viruses$keep_score_ind_vs2+viruses$keep_score_ind_ntv)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodmlydXNlcyRrZWVwX3Njb3JlX2luZF92czIsIHZpcnVzZXMka2VlcF9zY29yZV9pbmRfZHZmK3ZpcnVzZXMka2VlcF9zY29yZV9pbmRfdnMyK3ZpcnVzZXMka2VlcF9zY29yZV9pbmRfbnR2KVxuYGBgIn0= -->

```r
table(viruses$keep_score_ind_vs2, viruses$keep_score_ind_dvf+viruses$keep_score_ind_vs2+viruses$keep_score_ind_ntv)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

with 0.5 rules: 1464 sequences predicted viral that weren't viral by vs2 on its own
without them: 1234 sequences added


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodmlydXNlcyRrZWVwX3Njb3JlX2luZF92czIsIHZpcnVzZXMka2VlcF9zY29yZV9pbmRfZHZmK3ZpcnVzZXMka2VlcF9zY29yZV9pbmRfdnMyK3ZpcnVzZXMka2VlcF9zY29yZV9pbmRfbnR2LFxuICAgICAgdmlydXNlcyR0cnVlX3ZpcnVzKVxuYGBgIn0= -->

```r
table(viruses$keep_score_ind_vs2, viruses$keep_score_ind_dvf+viruses$keep_score_ind_vs2+viruses$keep_score_ind_ntv,
      viruses$true_virus)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

with 0.5 rules: however, of these, only 432 (30%) are truly viruses
without 0.5 rules: only 375 (30%) are truly viruses


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodmlydXNlcyRrZWVwX3Njb3JlX2luZF92czIsIHZpcnVzZXMka2VlcF9zY29yZV9pbmRfdmIrdmlydXNlcyRrZWVwX3Njb3JlX2luZF92czIrdmlydXNlcyRrZWVwX3Njb3JlX2luZF9udHYsXG4gICAgICB2aXJ1c2VzJHRydWVfdmlydXMpXG5gYGAifQ== -->

```r
table(viruses$keep_score_ind_vs2, viruses$keep_score_ind_vb+viruses$keep_score_ind_vs2+viruses$keep_score_ind_ntv,
      viruses$true_virus)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->

with 0.5 rules: however, of these, only 74 (7%) are truly viruses
without 0.5 rules: only 87 (4%) are truly viruses 



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucGFsX3B1cnBsZSA8LSBSQ29sb3JCcmV3ZXI6OmJyZXdlci5wYWwoNiwgXCJQdXJwbGVzXCIpXG5wYWxfb3JhbmdlIDwtIFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg0LCBcIk9yYW5nZXNcIilcblxuZ2dwbG90KHZpcnVzZXNfbG9uZ19zY29yZXMsIGFlcyh5PW4sIHg9cnVsZV9jb21iaW5hdGlvbixcbiAgICAgICAgICAgICAgICAgICBmaWxsPXZpcmFsX3Njb3JlKSkgK1xuICBnZW9tX2JhcihzdGF0PVwiaWRlbnRpdHlcIiwgY29sb3I9XCJibGFja1wiKSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICBjb29yZF9mbGlwKCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGxlZ2VuZC5wb3NpdGlvbiA9IFwiYm90dG9tXCIsXG4gICAgYXhpcy50ZXh0Lng9ZWxlbWVudF90ZXh0KHNpemU9OCwgYW5nbGUgPSA0NSlcbiAgKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKGMocmV2KHBhbF9vcmFuZ2UpLCBwYWxfcHVycGxlKSwgMC44KSkgK1xuICB4bGFiKFwiXCIpICtcbiAgeWxhYihcIk51bWJlciBvZiBTZXF1ZW5jZXNcIikgKyBcbiAgZmFjZXRfZ3JpZCh+c2VxdHlwZSwgc2NhbGVzID0gXCJmcmVlXCIpICtcbiAgc2NhbGVfeV9jb250aW51b3VzKGJyZWFrcyA9IHNjYWxlczo6cHJldHR5X2JyZWFrcyhuID0gMikpICtcbiAgc2NhbGVfeF9kaXNjcmV0ZShsYWJlbHM9cmV2KGMoXCJWaXJTb3J0ZXIyXCIsIFwiVmlyU29ydGVyXCIsIFwiVklCUkFOVFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwidHVuaW5nIHZpcmFsXCIsIFwidHVuaW5nIG5vdCB2aXJhbFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiRGVlcFZpckZpbmRlclwiKSkpXG5gYGAifQ== -->

```r
pal_purple <- RColorBrewer::brewer.pal(6, "Purples")
pal_orange <- RColorBrewer::brewer.pal(4, "Oranges")

ggplot(viruses_long_scores, aes(y=n, x=rule_combination,
                   fill=viral_score)) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    axis.text.x=element_text(size=8, angle = 45)
  ) +
  scale_fill_manual(name="",
                     values = alpha(c(rev(pal_orange), pal_purple), 0.8)) +
  xlab("") +
  ylab("Number of Sequences") + 
  facet_grid(~seqtype, scales = "free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  scale_x_discrete(labels=rev(c("VirSorter2", "VirSorter", "VIBRANT",
                            "tuning viral", "tuning not viral",
                            "DeepVirFinder")))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubnVtX2hpZ2hfdmlydXNlcyA8LSBtYXRyaXgoZGF0YT1OQSwgbnJvdz02MywgbmNvbD02MylcbnByb3BfaGlnaF92aXJ1c2VzIDwtIG1hdHJpeChkYXRhPU5BLCBucm93PTYzLCBuY29sPTYzKVxuXG5mb3IgKGkgaW4gMTo2Mykge1xuICBmb3IgKGogaW4gMTo2Mykge1xuICAgIG51bV9oaWdoX3ZpcnVzZXNbaSxqXSA8LSBzdW0odmlyYWxfc2NvcmVzWyxpXT49MSAmIHZpcmFsX3Njb3Jlc1ssal0+PTEpXG4gICAgXG4gICAgcHJvcF9oaWdoX3ZpcnVzZXNbaSxqXSA8LSBzdW0odmlyYWxfc2NvcmVzWyxpXT49MSAmIHZpcmFsX3Njb3Jlc1ssal0+PTEpXG4gICAgcHJvcF9oaWdoX3ZpcnVzZXNbaSxqXSA8LSBwcm9wX2hpZ2hfdmlydXNlc1tpLGpdL3N1bSh2aXJhbF9zY29yZXNbLGldPj0xIHwgdmlyYWxfc2NvcmVzWyxqXT49MSlcbiAgfVxufVxuXG5yb3duYW1lcyhudW1faGlnaF92aXJ1c2VzKSA8LSBjb21ib3NfbGlzdCR0b29sY29tYm9cbmNvbG5hbWVzKG51bV9oaWdoX3ZpcnVzZXMpIDwtIGNvbWJvc19saXN0JHRvb2xjb21ib1xuXG5yb3duYW1lcyhwcm9wX2hpZ2hfdmlydXNlcykgPC0gY29tYm9zX2xpc3QkdG9vbGNvbWJvXG5jb2xuYW1lcyhwcm9wX2hpZ2hfdmlydXNlcykgPC0gY29tYm9zX2xpc3QkdG9vbGNvbWJvXG5gYGAifQ== -->

```r
num_high_viruses <- matrix(data=NA, nrow=63, ncol=63)
prop_high_viruses <- matrix(data=NA, nrow=63, ncol=63)

for (i in 1:63) {
  for (j in 1:63) {
    num_high_viruses[i,j] <- sum(viral_scores[,i]>=1 & viral_scores[,j]>=1)
    
    prop_high_viruses[i,j] <- sum(viral_scores[,i]>=1 & viral_scores[,j]>=1)
    prop_high_viruses[i,j] <- prop_high_viruses[i,j]/sum(viral_scores[,i]>=1 | viral_scores[,j]>=1)
  }
}

rownames(num_high_viruses) <- combos_list$toolcombo
colnames(num_high_viruses) <- combos_list$toolcombo

rownames(prop_high_viruses) <- combos_list$toolcombo
colnames(prop_high_viruses) <- combos_list$toolcombo
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## visualizing a select set of tools

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxudGFibGUodmlydXNlcyRrZWVwX3Njb3JlX2luZF92czIsIHZpcnVzZXMka2VlcF9zY29yZV9pbmRfZHZmK3ZpcnVzZXMka2VlcF9zY29yZV9pbmRfdnMyK3ZpcnVzZXMka2VlcF9zY29yZV9pbmRfbnR2KVxuYGBgXG5gYGAifQ== -->

```r
```r
table(viruses$keep_score_ind_vs2, viruses$keep_score_ind_dvf+viruses$keep_score_ind_vs2+viruses$keep_score_ind_ntv)
```
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





# Comparing all of the tool methods - Fig 2A (also needed for 5)

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY29tYm9zX2xpc3QgPC0gZGF0YS5mcmFtZSh0b29sY29tYm89cmVwKDAsIDY0KSxcbmBgYCJ9 -->

```r
combos_list <- data.frame(toolcombo=rep(0, 64),
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiVGhlcmUgd2VyZSA0OCB3YXJuaW5ncyAodXNlIHdhcm5pbmdzKCkgdG8gc2VlIHRoZW0pXG4ifQ== -->

```
There were 48 warnings (use warnings() to see them)
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuICAgICAgICAgICAgICAgICAgICAgICAgICB0dW5lX25vdF92aXJhbD1yZXAoMCwgNjQpLFxuICAgICAgICAgICAgICAgICAgICAgICAgICBEVkY9cmVwKDAsIDY0KSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgdHVuZV92aXJhbD1yZXAoMCwgNjQpLFxuICAgICAgICAgICAgICAgICAgICAgICAgICBWSUJSQU5UPXJlcCgwLCA2NCksXG4gICAgICAgICAgICAgICAgICAgICAgICAgIFZTPXJlcCgwLCA2NCksXG4gICAgICAgICAgICAgICAgICAgICAgICAgIFZTMj1yZXAoMCwgNjQpKVxucCA8LSAxXG5cbmZvciAoaSBpbiBjKDAsMSkpe1xuICBmb3IgKGogaW4gYygwLDEpKXtcbiAgICBmb3IgKGsgaW4gYygwLDEpKXtcbiAgICAgIGZvciAobCBpbiBjKDAsMSkpe1xuICAgICAgICBmb3IgKG0gaW4gYygwLDEpKXtcbiAgICAgICAgICBmb3IgKG4gaW4gYygwLDEpKXtcbiAgICAgICAgICAgIGNvbWJvc19saXN0JHRvb2xjb21ib1twXSA8LSBwYXN0ZShpLGosayxsLG0sbilcbiAgICAgICAgICAgIGNvbWJvc19saXN0JHRvb2xjb21ibzJbcF0gPC0gcGFzdGUoaWYoaSl7XCJ0bnZcIn1lbHNle1wiMFwifSxpZihqKXtcIkRWRlwifWVsc2V7XCIwXCJ9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZihrKXtcInR2XCJ9ZWxzZXtcIjBcIn0saWYobCl7XCJWQlwifWVsc2V7XCIwXCJ9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZihtKXtcIlZTXCJ9ZWxzZXtcIjBcIn0saWYobil7XCJWUzJcIn1lbHNle1wiMFwifSlcbiAgICAgICAgICAgIGNvbWJvc19saXN0JHR1bmVfbm90X3ZpcmFsW3BdIDwtIGlcbiAgICAgICAgICAgIGNvbWJvc19saXN0JERWRltwXSA8LSBqXG4gICAgICAgICAgICBjb21ib3NfbGlzdCR0dW5lX3ZpcmFsW3BdIDwtIGtcbiAgICAgICAgICAgIGNvbWJvc19saXN0JFZJQlJBTlRbcF0gPC0gbFxuICAgICAgICAgICAgY29tYm9zX2xpc3QkVlNbcF0gPC0gbVxuICAgICAgICAgICAgY29tYm9zX2xpc3QkVlMyW3BdIDwtIG5cbiAgICAgICAgICAgIHAgPC0gcCsxXG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICB9XG59XG5cbmNvbWJvc19saXN0IDwtIGNvbWJvc19saXN0Wy0xLF1cblxudmlyYWxfc2NvcmVzIDwtIG1hdHJpeChkYXRhPTAsIG5yb3c9bnJvdyh2aXJ1c2VzKSwgbmNvbD1ucm93KGNvbWJvc19saXN0KSlcbm51bV92aXJ1c2VzIDwtIGRhdGEuZnJhbWUodG9vbGNvbWJvPXJlcCgwLCBucm93KGNvbWJvc19saXN0KSksXG4gICAgICAgICAgICAgICAgICAgICAgICAgIG51bV92aXJ1c2VzPXJlcCgwLCBucm93KGNvbWJvc19saXN0KSkpXG5cbmZvciAoaSBpbiAxOm5yb3coY29tYm9zX2xpc3QpKSB7XG4gIHZpcmFsX3Njb3Jlc1ssaV0gPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzLCBpbmNsdWRlX3ZpYnJhbnQgPSBjb21ib3NfbGlzdCRWSUJSQU5UW2ldLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlciA9IGNvbWJvc19saXN0JFZTW2ldLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBjb21ib3NfbGlzdCRWUzJbaV0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gY29tYm9zX2xpc3QkdHVuZV92aXJhbFtpXSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfbm90X3ZpcmFsID0gY29tYm9zX2xpc3QkdHVuZV9ub3RfdmlyYWxbaV0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IGNvbWJvc19saXN0JERWRltpXSlcbiAgXG4gIGlmIChtYXgodmlyYWxfc2NvcmVzWyxpXSk8PTApIHtcbiAgICBudW1fdmlydXNlcyRudW1fdmlydXNlc1tpXSA8LSAwXG4gIH1cbiAgZWxzZSB7XG4gICAgbnVtX3ZpcnVzZXMkbnVtX3ZpcnVzZXNbaV0gPC0gdGFibGUodmlyYWxfc2NvcmVzWyxpXT49MSlbWzJdXVxuICB9XG4gIG51bV92aXJ1c2VzJHRvb2xjb21ib1tpXSA8LSBjb21ib3NfbGlzdCR0b29sY29tYm9baV1cbiAgXG4gIG51bV92aXJ1c2VzJHRvb2xjb21ibzJbaV0gPC0gY29tYm9zX2xpc3QkdG9vbGNvbWJvMltpXVxufVxubnVtX3ZpcnVzZXMkbnVtdG9vbHMgPC0gc3RyX2NvdW50KG51bV92aXJ1c2VzJHRvb2xjb21ibywgXCIxXCIpXG5udW1fdmlydXNlcyA8LSBudW1fdmlydXNlc1tvcmRlcihudW1fdmlydXNlcyRudW1fdmlydXNlcywgZGVjcmVhc2luZz1GKSxdXG5udW1fdmlydXNlcyR0b29sY29tYm8gPC0gZmFjdG9yKG51bV92aXJ1c2VzJHRvb2xjb21ibywgbGV2ZWxzID0gdW5pcXVlKG51bV92aXJ1c2VzJHRvb2xjb21ibykpXG5udW1fdmlydXNlcyR0b29sY29tYm8yIDwtIGZhY3RvcihudW1fdmlydXNlcyR0b29sY29tYm8yLCBsZXZlbHMgPSB1bmlxdWUobnVtX3ZpcnVzZXMkdG9vbGNvbWJvMikpXG5udW1fdmlydXNlcyRudW10b29scyA8LSBhcy5mYWN0b3IobnVtX3ZpcnVzZXMkbnVtdG9vbHMpXG5cbnBhbCA8LSBnZ3RoZW1lczo6dGFibGVhdV9jb2xvcl9wYWwocGFsZXR0ZT1cIlRhYmxlYXUgMTBcIiwgdHlwZT1cInJlZ3VsYXJcIilcbnBhbDIgPC0gZ2d0aGVtZXM6OnRhYmxlYXVfY29sb3JfcGFsKHBhbGV0dGU9XCJUYWJsZWF1IDIwXCIsIHR5cGU9XCJyZWd1bGFyXCIpXG5cbnZpcnVzZXNfMSA8LSB2aXJ1c2VzW3ZpcnVzZXMkSW5kZXg9PVwiMVwiLF1cblxuYW5uX2NvbCA9IGRhdGEuZnJhbWUoXG4gICAgdG52ID0gY29tYm9zX2xpc3QkdHVuZV9ub3RfdmlyYWwsXG4gICAgdHYgPSBjb21ib3NfbGlzdCR0dW5lX3ZpcmFsLFxuICAgIGR2ZiA9IGNvbWJvc19saXN0JERWRixcbiAgICB2YiA9IGNvbWJvc19saXN0JFZJQlJBTlQsXG4gICAgdnMgPSBjb21ib3NfbGlzdCRWUyxcbiAgICB2czIgPSBjb21ib3NfbGlzdCRWUzJcbilcblxuYW5uX3JvdyA9IGRhdGEuZnJhbWUoXG4gICAgc2NvcmUgPSB2aXJ1c2VzXzEkc2VxdHlwZVxuKVxuXG5hbm5fY29sb3JzID0gbGlzdChcbiAgdG52ID0gYyhcImJsYWNrXCIsIHBhbCg2KVsxXSksXG4gIHR2ID0gYyhcImJsYWNrXCIsIHBhbCg2KVszXSksXG4gIGR2ZiA9IGMoXCJibGFja1wiLCBwYWwoNilbMl0pLFxuICB2YiA9IGMoXCJibGFja1wiLCBwYWwoNilbNF0pLFxuICB2cyA9IGMoXCJibGFja1wiLCBwYWwoNilbNV0pLFxuICB2czIgPSBjKFwiYmxhY2tcIiwgcGFsKDYpWzZdKSxcbiAgc2NvcmUgPSBjKFwiYXJjaGFlYVwiPXBhbDIoMjApWzEzXSwgXCJiYWN0ZXJpYVwiPXBhbDIoMjApWzldLCBcImZ1bmdpXCI9cGFsMigyMClbMjBdLFxuICAgICAgICAgICAgICBcInBsYXNtaWRcIj1wYWwyKDIwKVsxOV0sIFwicHJvdGlzdFwiPXBhbDIoMjApWzVdLCBcInZpcnVzXCI9cGFsMigyMClbMTddKVxuKVxuXG5yb3duYW1lcyhhbm5fY29sKSA8LSBjb21ib3NfbGlzdCR0b29sY29tYm9cbnJvd25hbWVzKGFubl9yb3cpIDwtIHZpcnVzZXNfMSRjb250aWdcblxudmlyYWxfc2NvcmVzXzEgPC0gdmlyYWxfc2NvcmVzW3ZpcnVzZXMkSW5kZXg9PVwiMVwiLF1cblxuY29sbmFtZXModmlyYWxfc2NvcmVzXzEpIDwtIGNvbWJvc19saXN0JHRvb2xjb21ib1xucm93bmFtZXModmlyYWxfc2NvcmVzXzEpIDwtIHZpcnVzZXNfMSRjb250aWdcblxudmlyYWxfc2NvcmVzXzEgPC0gZmxvb3IodmlyYWxfc2NvcmVzXzEpXG5cbnNlcV9zY29yZXMgPC0gcm93U3Vtcyh2aXJhbF9zY29yZXNfMSlcblxuc2VxX3R5cGUgPC0gdGFibGUodmlydXNlc18xJHNlcXR5cGUpW2MoNiwyLDEsNCw1LDMpXVxuc2VxX3R5cGUgPC0gY3Vtc3VtKHNlcV90eXBlKVxuXG5zdWJfc2VxX3Njb3JlcyA8LSBzZXFfc2NvcmVzWzE6c2VxX3R5cGVbMV1dXG5vcmRlcl9zZXFfc2NvcmVzX3ZpcnVzIDwtIG5hbWVzKHN1Yl9zZXFfc2NvcmVzKVtvcmRlcihzdWJfc2VxX3Njb3JlcyldXG5cbnN1Yl9zZXFfc2NvcmVzIDwtIHNlcV9zY29yZXNbc2VxX3R5cGVbMV06c2VxX3R5cGVbMl1dXG5vcmRlcl9zZXFfc2NvcmVzX2JhY3RlcmlhIDwtIG5hbWVzKHN1Yl9zZXFfc2NvcmVzKVtvcmRlcihzdWJfc2VxX3Njb3JlcyldXG5cbnN1Yl9zZXFfc2NvcmVzIDwtIHNlcV9zY29yZXNbc2VxX3R5cGVbMl06c2VxX3R5cGVbM11dXG5vcmRlcl9zZXFfc2NvcmVzX2FyY2hhZWEgPC0gbmFtZXMoc3ViX3NlcV9zY29yZXMpW29yZGVyKHN1Yl9zZXFfc2NvcmVzKV1cblxuc3ViX3NlcV9zY29yZXMgPC0gc2VxX3Njb3Jlc1tzZXFfdHlwZVszXTpzZXFfdHlwZVs0XV1cbm9yZGVyX3NlcV9zY29yZXNfcGxhc21pZCA8LSBuYW1lcyhzdWJfc2VxX3Njb3Jlcylbb3JkZXIoc3ViX3NlcV9zY29yZXMpXVxuXG5zdWJfc2VxX3Njb3JlcyA8LSBzZXFfc2NvcmVzW3NlcV90eXBlWzRdOnNlcV90eXBlWzVdXVxub3JkZXJfc2VxX3Njb3Jlc19wcm90aXN0IDwtIG5hbWVzKHN1Yl9zZXFfc2NvcmVzKVtvcmRlcihzdWJfc2VxX3Njb3JlcyldXG5cbnN1Yl9zZXFfc2NvcmVzIDwtIHNlcV9zY29yZXNbc2VxX3R5cGVbNV06c2VxX3R5cGVbNl1dXG5vcmRlcl9zZXFfc2NvcmVzX2Z1bmdpIDwtIG5hbWVzKHN1Yl9zZXFfc2NvcmVzKVtvcmRlcihzdWJfc2VxX3Njb3JlcyldXG5cbnZpcmFsX3Njb3Jlc18xIDwtIHZpcmFsX3Njb3Jlc18xW2Mob3JkZXJfc2VxX3Njb3Jlc192aXJ1cywgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG9yZGVyX3NlcV9zY29yZXNfYmFjdGVyaWEsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG9yZGVyX3NlcV9zY29yZXNfYXJjaGFlYSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgb3JkZXJfc2VxX3Njb3Jlc19wbGFzbWlkLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBvcmRlcl9zZXFfc2NvcmVzX3Byb3Rpc3QsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG9yZGVyX3NlcV9zY29yZXNfZnVuZ2kpLF1cblxucGhlYXRtYXA6OnBoZWF0bWFwKHQodmlyYWxfc2NvcmVzXzEpLFxuICAgICAgICAgICAgICAgICAgYW5ub3RhdGlvbl9jb2wgPSBhbm5fcm93LFxuICAgICAgICAgICAgICAgICAgYW5ub3RhdGlvbl9yb3cgPSBhbm5fY29sLFxuICAgICAgICAgICAgICAgICAgYW5ub3RhdGlvbl9jb2xvcnMgPSBhbm5fY29sb3JzLFxuICAgICAgICAgICAgICAgICAgI2JvcmRlcl9jb2xvciA9IFwiYmxhY2tcIixcbiAgICAgICAgICAgICAgICAgIGNvbG9yID0gYyhyZXYoUkNvbG9yQnJld2VyOjpicmV3ZXIucGFsKDksXCJPcmFuZ2VzXCIpWzQ6N10pLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIFJDb2xvckJyZXdlcjo6YnJld2VyLnBhbCg5LFwiUHVycGxlc1wiKVs0OjldIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICksXG4gICAgICAgICAgICAgICAgICBzaG93X3Jvd25hbWVzID0gRixcbiAgICAgICAgICAgICAgICAgIHNob3dfY29sbmFtZXMgPSBGLFxuICAgICAgICAgICAgICAgICAgYW5ub3RhdGlvbl9sZWdlbmQgPSBGLFxuICAgICAgICAgICAgICAgICAgY2x1c3Rlcl9yb3dzID0gRixcbiAgICAgICAgICAgICAgICAgIGNsdXN0ZXJfY29scyA9IEYsXG4gICAgICAgICAgICAgICAgICBjZWxsaGVpZ2h0PTMsXG4gICAgICAgICAgICAgICAgICBnYXBzX2NvbCA9IGMoc2VxX3R5cGVbMV0sIHNlcV90eXBlWzJdLCBzZXFfdHlwZVszXSwgc2VxX3R5cGVbNF0sIHNlcV90eXBlWzVdKSxcbiAgICAgICAgICAgICAgICAgIGxlZ2VuZF9icmVha3M9YygtMzo1KSxcbiAgICAgICAgICAgICAgICAgIGxlZ2VuZF9sYWJlbHM9YygtMzo1KSxcbiAgICAgICAgICAgICAgICAgIGZpbGVuYW1lID0gXCIuLi9JbnRlcm1lZGlhcnlGaWxlcy9waGVhdG1hcF92aXJ1c19zZXRfMV92aXJhbF9zY29yZXMucG5nXCJcbiAgICAgICAgICAgICAgICAgIClcbmBgYCJ9 -->

```r
                          tune_not_viral=rep(0, 64),
                          DVF=rep(0, 64),
                          tune_viral=rep(0, 64),
                          VIBRANT=rep(0, 64),
                          VS=rep(0, 64),
                          VS2=rep(0, 64))
p <- 1

for (i in c(0,1)){
  for (j in c(0,1)){
    for (k in c(0,1)){
      for (l in c(0,1)){
        for (m in c(0,1)){
          for (n in c(0,1)){
            combos_list$toolcombo[p] <- paste(i,j,k,l,m,n)
            combos_list$toolcombo2[p] <- paste(if(i){"tnv"}else{"0"},if(j){"DVF"}else{"0"},
                                               if(k){"tv"}else{"0"},if(l){"VB"}else{"0"},
                                               if(m){"VS"}else{"0"},if(n){"VS2"}else{"0"})
            combos_list$tune_not_viral[p] <- i
            combos_list$DVF[p] <- j
            combos_list$tune_viral[p] <- k
            combos_list$VIBRANT[p] <- l
            combos_list$VS[p] <- m
            combos_list$VS2[p] <- n
            p <- p+1
          }
        }
      }
    }
  }
}

combos_list <- combos_list[-1,]

viral_scores <- matrix(data=0, nrow=nrow(viruses), ncol=nrow(combos_list))
num_viruses <- data.frame(toolcombo=rep(0, nrow(combos_list)),
                          num_viruses=rep(0, nrow(combos_list)))

for (i in 1:nrow(combos_list)) {
  viral_scores[,i] <- getting_viral_set_1(viruses, include_vibrant = combos_list$VIBRANT[i],
                                            include_virsorter = combos_list$VS[i],
                                            include_virsorter2 = combos_list$VS2[i],
                                            include_tuning_viral = combos_list$tune_viral[i],
                                            include_tuning_not_viral = combos_list$tune_not_viral[i],
                                            include_deepvirfinder = combos_list$DVF[i])
  
  if (max(viral_scores[,i])<=0) {
    num_viruses$num_viruses[i] <- 0
  }
  else {
    num_viruses$num_viruses[i] <- table(viral_scores[,i]>=1)[[2]]
  }
  num_viruses$toolcombo[i] <- combos_list$toolcombo[i]
  
  num_viruses$toolcombo2[i] <- combos_list$toolcombo2[i]
}
num_viruses$numtools <- str_count(num_viruses$toolcombo, "1")
num_viruses <- num_viruses[order(num_viruses$num_viruses, decreasing=F),]
num_viruses$toolcombo <- factor(num_viruses$toolcombo, levels = unique(num_viruses$toolcombo))
num_viruses$toolcombo2 <- factor(num_viruses$toolcombo2, levels = unique(num_viruses$toolcombo2))
num_viruses$numtools <- as.factor(num_viruses$numtools)

pal <- ggthemes::tableau_color_pal(palette="Tableau 10", type="regular")
pal2 <- ggthemes::tableau_color_pal(palette="Tableau 20", type="regular")

viruses_1 <- viruses[viruses$Index=="1",]

ann_col = data.frame(
    tnv = combos_list$tune_not_viral,
    tv = combos_list$tune_viral,
    dvf = combos_list$DVF,
    vb = combos_list$VIBRANT,
    vs = combos_list$VS,
    vs2 = combos_list$VS2
)

ann_row = data.frame(
    score = viruses_1$seqtype
)

ann_colors = list(
  tnv = c("black", pal(6)[1]),
  tv = c("black", pal(6)[3]),
  dvf = c("black", pal(6)[2]),
  vb = c("black", pal(6)[4]),
  vs = c("black", pal(6)[5]),
  vs2 = c("black", pal(6)[6]),
  score = c("archaea"=pal2(20)[13], "bacteria"=pal2(20)[9], "fungi"=pal2(20)[20],
              "plasmid"=pal2(20)[19], "protist"=pal2(20)[5], "virus"=pal2(20)[17])
)

rownames(ann_col) <- combos_list$toolcombo
rownames(ann_row) <- viruses_1$contig

viral_scores_1 <- viral_scores[viruses$Index=="1",]

colnames(viral_scores_1) <- combos_list$toolcombo
rownames(viral_scores_1) <- viruses_1$contig

viral_scores_1 <- floor(viral_scores_1)

seq_scores <- rowSums(viral_scores_1)

seq_type <- table(viruses_1$seqtype)[c(6,2,1,4,5,3)]
seq_type <- cumsum(seq_type)

sub_seq_scores <- seq_scores[1:seq_type[1]]
order_seq_scores_virus <- names(sub_seq_scores)[order(sub_seq_scores)]

sub_seq_scores <- seq_scores[seq_type[1]:seq_type[2]]
order_seq_scores_bacteria <- names(sub_seq_scores)[order(sub_seq_scores)]

sub_seq_scores <- seq_scores[seq_type[2]:seq_type[3]]
order_seq_scores_archaea <- names(sub_seq_scores)[order(sub_seq_scores)]

sub_seq_scores <- seq_scores[seq_type[3]:seq_type[4]]
order_seq_scores_plasmid <- names(sub_seq_scores)[order(sub_seq_scores)]

sub_seq_scores <- seq_scores[seq_type[4]:seq_type[5]]
order_seq_scores_protist <- names(sub_seq_scores)[order(sub_seq_scores)]

sub_seq_scores <- seq_scores[seq_type[5]:seq_type[6]]
order_seq_scores_fungi <- names(sub_seq_scores)[order(sub_seq_scores)]

viral_scores_1 <- viral_scores_1[c(order_seq_scores_virus, 
                                   order_seq_scores_bacteria,
                                   order_seq_scores_archaea,
                                   order_seq_scores_plasmid,
                                   order_seq_scores_protist,
                                   order_seq_scores_fungi),]

pheatmap::pheatmap(t(viral_scores_1),
                  annotation_col = ann_row,
                  annotation_row = ann_col,
                  annotation_colors = ann_colors,
                  #border_color = "black",
                  color = c(rev(RColorBrewer::brewer.pal(9,"Oranges")[4:7]),
                            RColorBrewer::brewer.pal(9,"Purples")[4:9] 
                            ),
                  show_rownames = F,
                  show_colnames = F,
                  annotation_legend = F,
                  cluster_rows = F,
                  cluster_cols = F,
                  cellheight=3,
                  gaps_col = c(seq_type[1], seq_type[2], seq_type[3], seq_type[4], seq_type[5]),
                  legend_breaks=c(-3:5),
                  legend_labels=c(-3:5),
                  filename = "../IntermediaryFiles/pheatmap_virus_set_1_viral_scores.png"
                  )
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


# comparing all rules against each other - Figure 5

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubnVtX2hpZ2hfdmlydXNlcyA8LSBtYXRyaXgoZGF0YT1OQSwgbnJvdz02MywgbmNvbD02MylcbmBgYCJ9 -->

```r
num_high_viruses <- matrix(data=NA, nrow=63, ncol=63)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiVGhlcmUgd2VyZSAxNiB3YXJuaW5ncyAodXNlIHdhcm5pbmdzKCkgdG8gc2VlIHRoZW0pXG4ifQ== -->

```
There were 16 warnings (use warnings() to see them)
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucHJvcF9oaWdoX3ZpcnVzZXMgPC0gbWF0cml4KGRhdGE9TkEsIG5yb3c9NjMsIG5jb2w9NjMpXG5cbmZvciAoaSBpbiAxOjYzKSB7XG4gIGZvciAoaiBpbiAxOjYzKSB7XG4gICAgbnVtX2hpZ2hfdmlydXNlc1tpLGpdIDwtIHN1bSh2aXJhbF9zY29yZXNbLGldPj0xICYgdmlyYWxfc2NvcmVzWyxqXT49MSlcbiAgICBcbiAgICBwcm9wX2hpZ2hfdmlydXNlc1tpLGpdIDwtIHN1bSh2aXJhbF9zY29yZXNbLGldPj0xICYgdmlyYWxfc2NvcmVzWyxqXT49MSlcbiAgICBwcm9wX2hpZ2hfdmlydXNlc1tpLGpdIDwtIHByb3BfaGlnaF92aXJ1c2VzW2ksal0vc3VtKHZpcmFsX3Njb3Jlc1ssaV0+PTEgfCB2aXJhbF9zY29yZXNbLGpdPj0xKVxuICB9XG59XG5cbnJvd25hbWVzKG51bV9oaWdoX3ZpcnVzZXMpIDwtIGNvbWJvc19saXN0JHRvb2xjb21ib1xuY29sbmFtZXMobnVtX2hpZ2hfdmlydXNlcykgPC0gY29tYm9zX2xpc3QkdG9vbGNvbWJvXG5cbnJvd25hbWVzKHByb3BfaGlnaF92aXJ1c2VzKSA8LSBjb21ib3NfbGlzdCR0b29sY29tYm9cbmNvbG5hbWVzKHByb3BfaGlnaF92aXJ1c2VzKSA8LSBjb21ib3NfbGlzdCR0b29sY29tYm9cbmBgYCJ9 -->

```r
prop_high_viruses <- matrix(data=NA, nrow=63, ncol=63)

for (i in 1:63) {
  for (j in 1:63) {
    num_high_viruses[i,j] <- sum(viral_scores[,i]>=1 & viral_scores[,j]>=1)
    
    prop_high_viruses[i,j] <- sum(viral_scores[,i]>=1 & viral_scores[,j]>=1)
    prop_high_viruses[i,j] <- prop_high_viruses[i,j]/sum(viral_scores[,i]>=1 | viral_scores[,j]>=1)
  }
}

rownames(num_high_viruses) <- combos_list$toolcombo
colnames(num_high_viruses) <- combos_list$toolcombo

rownames(prop_high_viruses) <- combos_list$toolcombo
colnames(prop_high_viruses) <- combos_list$toolcombo
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucGhlYXRtYXA6OnBoZWF0bWFwKHByb3BfaGlnaF92aXJ1c2VzLFxuXG5gYGAifQ== -->

```r
pheatmap::pheatmap(prop_high_viruses,

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiVGhlcmUgd2VyZSAxNiB3YXJuaW5ncyAodXNlIHdhcm5pbmdzKCkgdG8gc2VlIHRoZW0pXG4ifQ== -->

```
There were 16 warnings (use warnings() to see them)
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuICAgICAgICAgICAgICAgICAgYW5ub3RhdGlvbl9yb3cgPSBhbm5fY29sLFxuICAgICAgICAgICAgICAgICAgYW5ub3RhdGlvbl9jb2xvcnMgPSBhbm5fY29sb3JzLFxuICAgICAgICAgICAgICAgICAgYm9yZGVyX2NvbG9yID0gXCJibGFja1wiLFxuICAgICAgICAgICAgICAgICAgY29sb3IgPSB2aXJpZGlzKDEwMCksXG4gICAgICAgICAgICAgICAgICBjdXRyZWVfcm93cyA9IDMsXG4gICAgICAgICAgICAgICAgICBjdXRyZWVfY29scyA9IDMsXG4gICAgICAgICAgICAgICAgICBzaG93X3Jvd25hbWVzID0gRixcbiAgICAgICAgICAgICAgICAgIHNob3dfY29sbmFtZXMgPSBGLFxuICAgICAgICAgICAgICAgICAgZmlsZW5hbWUgPSBcIi4uL0ludGVybWVkaWFyeUZpbGVzL3BoZWF0bWFwX3Byb3BfdmlydXNlc19pbl9jb21tb24ucG5nXCIsXG4gICAgICAgICAgICAgICAgICBhbm5vdGF0aW9uX2xlZ2VuZCA9IEYpXG5gYGAifQ== -->

```r
                  annotation_row = ann_col,
                  annotation_colors = ann_colors,
                  border_color = "black",
                  color = viridis(100),
                  cutree_rows = 3,
                  cutree_cols = 3,
                  show_rownames = F,
                  show_colnames = F,
                  filename = "../IntermediaryFiles/pheatmap_prop_viruses_in_common.png",
                  annotation_legend = F)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





# numbers in common between two tool combos

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY29tYm9zX2xpc3RfdHdvIDwtIGNvbWJvc19saXN0XG5jb21ib3NfbGlzdF90d28kbnVtdG9vbHMgPC0gc3RyX2NvdW50KGNvbWJvc19saXN0X3R3byR0b29sY29tYm8sIFwiMVwiKVxuY29tYm9zX2xpc3RfdHdvIDwtIGNvbWJvc19saXN0X3R3b1tjb21ib3NfbGlzdF90d28kbnVtdG9vbHM9PTEsXVxuXG50d29fdG9vbF9wYWlycyA8LSBudW1faGlnaF92aXJ1c2VzW3Jvd25hbWVzKG51bV9oaWdoX3ZpcnVzZXMpICVpbiUgY29tYm9zX2xpc3RfdHdvJHRvb2xjb21ibyxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sbmFtZXMobnVtX2hpZ2hfdmlydXNlcykgJWluJSBjb21ib3NfbGlzdF90d28kdG9vbGNvbWJvXVxuXG5jb2xuYW1lcyh0d29fdG9vbF9wYWlycykgPC0gYyhcIlZTMlwiLCBcIlZTXCIsIFwiVkJcIiwgXCJ0dlwiLCBcIkRWRlwiLCBcInRudlwiKVxucm93bmFtZXModHdvX3Rvb2xfcGFpcnMpIDwtIGMoXCJWUzJcIiwgXCJWU1wiLCBcIlZCXCIsIFwidHZcIiwgXCJEVkZcIiwgXCJ0bnZcIilcbmBgYCJ9 -->

```r
combos_list_two <- combos_list
combos_list_two$numtools <- str_count(combos_list_two$toolcombo, "1")
combos_list_two <- combos_list_two[combos_list_two$numtools==1,]

two_tool_pairs <- num_high_viruses[rownames(num_high_viruses) %in% combos_list_two$toolcombo,
                                   colnames(num_high_viruses) %in% combos_list_two$toolcombo]

colnames(two_tool_pairs) <- c("VS2", "VS", "VB", "tv", "DVF", "tnv")
rownames(two_tool_pairs) <- c("VS2", "VS", "VB", "tv", "DVF", "tnv")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlyYWxfc2NvcmVzX3R3byA8LSB2aXJhbF9zY29yZXNbLGNvbG5hbWVzKG51bV9oaWdoX3ZpcnVzZXMpICVpbiUgY29tYm9zX2xpc3RfdHdvJHRvb2xjb21ib11cblxubiA8LSA2XG52aXJhbF9zY29yZXNfdHdvX2NvbWIgPC0gbWF0cml4KGRhdGE9TkEsIG5yb3c9biwgbmNvbD1uKVxuXG5mb3IgKGkgaW4gMTpuKSB7XG4gIGZvciAoaiBpbiAxOm4pIHtcbiAgICBpZiAoaT09aikge1xuICAgICAgdmlyYWxfc2NvcmVzX3R3b19jb21iW2ksal0gPC0gc3VtKHZpcmFsX3Njb3Jlc190d29bLGldPj0xKVxuICAgIH1cbiAgICBlbHNlIGlmIChpPGopIHtcbiAgICAgIHZpcmFsX3Njb3Jlc190d29fY29tYltpLGpdIDwtIHN1bSgodmlyYWxfc2NvcmVzX3R3b1ssaV0rdmlyYWxfc2NvcmVzX3R3b1ssal0pPj0xKVxuICAgIH1cbiAgICBlbHNlIGlmIChpPmopIHtcbiAgICAgIHZpcmFsX3Njb3Jlc190d29fY29tYltpLGpdIDwtIHBhc3RlKHJvdW5kKHN1bSgodmlyYWxfc2NvcmVzX3R3b1ssaV0rdmlyYWxfc2NvcmVzX3R3b1ssal0pPj0xKS9zdW0odmlyYWxfc2NvcmVzX3R3b1ssaV0+PTEgfCB2aXJhbF9zY29yZXNfdHdvWyxqXT49MSkqMTAwLCBkaWdpdHM9MCksIFwiJVwiLCBzZXA9XCJcIilcbiAgICB9XG4gIH1cbn1cblxuY29sbmFtZXModmlyYWxfc2NvcmVzX3R3b19jb21iKSA8LSBjKFwiVlMyXCIsIFwiVlNcIiwgXCJWQlwiLCBcInR2XCIsIFwiRFZGXCIsIFwidG52XCIpXG5yb3duYW1lcyh2aXJhbF9zY29yZXNfdHdvX2NvbWIpIDwtIGMoXCJWUzJcIiwgXCJWU1wiLCBcIlZCXCIsIFwidHZcIiwgXCJEVkZcIiwgXCJ0bnZcIilcbmBgYCJ9 -->

```r
viral_scores_two <- viral_scores[,colnames(num_high_viruses) %in% combos_list_two$toolcombo]

n <- 6
viral_scores_two_comb <- matrix(data=NA, nrow=n, ncol=n)

for (i in 1:n) {
  for (j in 1:n) {
    if (i==j) {
      viral_scores_two_comb[i,j] <- sum(viral_scores_two[,i]>=1)
    }
    else if (i<j) {
      viral_scores_two_comb[i,j] <- sum((viral_scores_two[,i]+viral_scores_two[,j])>=1)
    }
    else if (i>j) {
      viral_scores_two_comb[i,j] <- paste(round(sum((viral_scores_two[,i]+viral_scores_two[,j])>=1)/sum(viral_scores_two[,i]>=1 | viral_scores_two[,j]>=1)*100, digits=0), "%", sep="")
    }
  }
}

colnames(viral_scores_two_comb) <- c("VS2", "VS", "VB", "tv", "DVF", "tnv")
rownames(viral_scores_two_comb) <- c("VS2", "VS", "VB", "tv", "DVF", "tnv")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlyYWxfc2NvcmVzX3R3b19jb21iXG5gYGAifQ== -->

```r
viral_scores_two_comb
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


comparing VS2 and VIBRANT viruses - Figure S15

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudnMyX3ZiIDwtIGRhdGEuZnJhbWUodnMyPXZpcmFsX3Njb3Jlc190d29bLDFdLFxuICAgICAgICAgICAgICAgICAgICAgdmI9dmlyYWxfc2NvcmVzX3R3b1ssM10sXG4gICAgICAgICAgICAgICAgICAgICBjb21iX3Njb3JlPXZpcmFsX3Njb3Jlc190d29bLDFdK3ZpcmFsX3Njb3Jlc190d29bLDNdLFxuICAgICAgICAgICAgICAgICAgICAgc2VxdHlwZT12aXJ1c2VzJHNlcXR5cGUsXG4gICAgICAgICAgICAgICAgICAgICBtYXhfc2NvcmVfZ3JvdXA9dmlydXNlcyRtYXhfc2NvcmVfZ3JvdXApXG5cbnZzMl92YiRjb21iW3ZzMl92YiRjb21iX3Njb3JlPDFdIDwtIFwibmVpdGhlclwiXG52czJfdmIkY29tYlt2czJfdmIkdnMyPT0wLjUgJiB2czJfdmIkdmI9PTAuNV0gPC0gXCJib3RoX2xvd1wiXG52czJfdmIkY29tYlt2czJfdmIkdnMyPj0xICYgdnMyX3ZiJHZiPT0wLjVdIDwtIFwidnMyX2hpZ2hfdmJfbG93XCJcbnZzMl92YiRjb21iW3ZzMl92YiR2czI+PTAuNSAmIHZzMl92YiR2Yj09MV0gPC0gXCJ2czJfbG93X3ZiX2hpZ2hcIlxudnMyX3ZiJGNvbWJbdnMyX3ZiJHZzMj49MSAmIHZzMl92YiR2Yj09MF0gPC0gXCJ2czJfb25seVwiXG52czJfdmIkY29tYlt2czJfdmIkdnMyPT0wICYgdnMyX3ZiJHZiPj0xXSA8LSBcInZiX29ubHlcIlxudnMyX3ZiJGNvbWJbdnMyX3ZiJHZzMj49MSAmIHZzMl92YiR2Yj49MV0gPC0gXCJib3RoX2hpZ2hcIlxuXG4jdnMyX3ZiJHRydWVfdmlydXMgPC0gXCJub3RcIlxuI3ZzMl92YiR0cnVlX3ZpcnVzW3ZzMl92YiRzZXF0eXBlPT1cInZpcnVzXCJdIDwtIFwidmlydXNcIlxuXG52czJfdmJfc2VxdHlwZSA8LSB2czJfdmIgJT4lIFxuICBncm91cF9ieShjb21iLCBjb21iX3Njb3JlLCBzZXF0eXBlKSAlPiVcbiAgc3VtbWFyaXNlKG4gPSBuKCkpXG5gYGAifQ== -->

```r
vs2_vb <- data.frame(vs2=viral_scores_two[,1],
                     vb=viral_scores_two[,3],
                     comb_score=viral_scores_two[,1]+viral_scores_two[,3],
                     seqtype=viruses$seqtype,
                     max_score_group=viruses$max_score_group)

vs2_vb$comb[vs2_vb$comb_score<1] <- "neither"
vs2_vb$comb[vs2_vb$vs2==0.5 & vs2_vb$vb==0.5] <- "both_low"
vs2_vb$comb[vs2_vb$vs2>=1 & vs2_vb$vb==0.5] <- "vs2_high_vb_low"
vs2_vb$comb[vs2_vb$vs2>=0.5 & vs2_vb$vb==1] <- "vs2_low_vb_high"
vs2_vb$comb[vs2_vb$vs2>=1 & vs2_vb$vb==0] <- "vs2_only"
vs2_vb$comb[vs2_vb$vs2==0 & vs2_vb$vb>=1] <- "vb_only"
vs2_vb$comb[vs2_vb$vs2>=1 & vs2_vb$vb>=1] <- "both_high"

#vs2_vb$true_virus <- "not"
#vs2_vb$true_virus[vs2_vb$seqtype=="virus"] <- "virus"

vs2_vb_seqtype <- vs2_vb %>% 
  group_by(comb, comb_score, seqtype) %>%
  summarise(n = n())
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUodnMyX3ZiJGNvbWIpXG5gYGAifQ== -->

```r
table(vs2_vb$comb)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KHZzMl92Yl9zZXF0eXBlLCBhZXMoeT1uLCB4PWNvbWIsXG4gICAgICAgICAgICAgICAgICAgZmlsbD1hcy5mYWN0b3IoY29tYl9zY29yZSkpKSArXG4gIGdlb21fYmFyKHN0YXQ9XCJpZGVudGl0eVwiLCBjb2xvcj1cImJsYWNrXCIpICtcbiAgdGhlbWVfbGlnaHQoKSArXG4gIGNvb3JkX2ZsaXAoKSArXG4gIHRoZW1lKFxuICAgIHBhbmVsLmdyaWQubWFqb3IueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIlxuICApICtcbiAgc2NhbGVfZmlsbF9icmV3ZXIocGFsZXR0ZSA9IFwiUHVPclwiLCApICtcbiAgeGxhYihcIlwiKSArXG4gIHlsYWIoXCJOdW1iZXIgb2YgU2VxdWVuY2VzXCIpICsgXG4gIGZhY2V0X2dyaWQofnNlcXR5cGUsIHNjYWxlcyA9IFwiZnJlZVwiKVxuYGBgIn0= -->

```r
ggplot(vs2_vb_seqtype, aes(y=n, x=comb,
                   fill=as.factor(comb_score))) +
  geom_bar(stat="identity", color="black") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_brewer(palette = "PuOr", ) +
  xlab("") +
  ylab("Number of Sequences") + 
  facet_grid(~seqtype, scales = "free")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


# overlap in proportion of contigs for tool sets - Figure S16

get number of tools and number of rules for each combination

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYW5uX2NvbCRudW1ydWxlcyA8LSByb3dTdW1zKGFubl9jb2xbLDE6Nl0pXG5hbm5fY29sJG51bXRvb2xzIDwtIGFubl9jb2wkbnVtcnVsZXNcblxuZm9yIChpIGluIDE6bnJvdyhhbm5fY29sKSkge1xuICBpZiAoYW5uX2NvbCR0dltpXT09MSkge1xuICAgIGFubl9jb2wkbnVtdG9vbHNbaV0gPC0gYW5uX2NvbCRudW10b29sc1tpXSArIDNcbiAgfVxuICBcbiAgaWYgKGFubl9jb2wkdG52W2ldPT0xKSB7XG4gICAgYW5uX2NvbCRudW10b29sc1tpXSA8LSBhbm5fY29sJG51bXRvb2xzW2ldICsgNVxuICB9XG59XG5cbmFubl9jb2wkbnVtdG9vbHNbYW5uX2NvbCRudW10b29scz42XSA8LSA2XG5cbmFubl9jb2wkdG9vbCA8LSByb3duYW1lcyhhbm5fY29sKVxuYGBgIn0= -->

```r
ann_col$numrules <- rowSums(ann_col[,1:6])
ann_col$numtools <- ann_col$numrules

for (i in 1:nrow(ann_col)) {
  if (ann_col$tv[i]==1) {
    ann_col$numtools[i] <- ann_col$numtools[i] + 3
  }
  
  if (ann_col$tnv[i]==1) {
    ann_col$numtools[i] <- ann_col$numtools[i] + 5
  }
}

ann_col$numtools[ann_col$numtools>6] <- 6

ann_col$tool <- rownames(ann_col)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucHJvcF9oaWdoX3ZpcnVzZXNfbG9uZyA8LSBhcy5kYXRhLmZyYW1lKHByb3BfaGlnaF92aXJ1c2VzKVxucHJvcF9oaWdoX3ZpcnVzZXNfbG9uZyR0b29sIDwtIGNvbWJvc19saXN0JHRvb2xjb21ib1xuXG5wcm9wX2hpZ2hfdmlydXNlc19sb25nIDwtIGFzLmRhdGEuZnJhbWUocHJvcF9oaWdoX3ZpcnVzZXNfbG9uZykgJT4lIFxuICBwaXZvdF9sb25nZXIoIXRvb2wsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XCJwcm9wb3J0aW9uXCIsXG4gICAgICAgICAgICAgICBuYW1lc190bz1cInRvb2wyXCIpICU+JVxuICBpbm5lcl9qb2luKGFubl9jb2wsIGJ5PWMoXCJ0b29sXCIpKVxuXG5cbmBgYCJ9 -->

```r
prop_high_viruses_long <- as.data.frame(prop_high_viruses)
prop_high_viruses_long$tool <- combos_list$toolcombo

prop_high_viruses_long <- as.data.frame(prop_high_viruses_long) %>% 
  pivot_longer(!tool,
               values_to="proportion",
               names_to="tool2") %>%
  inner_join(ann_col, by=c("tool"))

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuaGVhZChwcm9wX2hpZ2hfdmlydXNlc19sb25nKVxuYGBgIn0= -->

```r
head(prop_high_viruses_long)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KHByb3BfaGlnaF92aXJ1c2VzX2xvbmdbcHJvcF9oaWdoX3ZpcnVzZXNfbG9uZyR2czI9PTEsXSwgYWVzKHg9YXMuZmFjdG9yKG51bXJ1bGVzKSwgeT1wcm9wb3J0aW9uLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb2xvcj1hcy5mYWN0b3IobnVtcnVsZXMpLCBmaWxsPWFzLmZhY3RvcihudW1ydWxlcykpKSArXG4gIGdlb21fYm94cGxvdCgpICtcbiAgI2dlb21fcG9pbnQoKSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCwgYW5nbGUgPSA5MCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTYpLFxuICApICtcbiAgeGxhYihcIk51bWJlciBvZiBSdWxlc1wiKSArXG4gIHlsYWIoXCJQcm9wb3J0aW9uIG9mIENvbnRpZ3MgaW4gQ29tbW9uXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAwLjMpKSArXG4gIHNjYWxlX2NvbG9yX21hbnVhbChuYW1lPVwiXCIsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShyZXYodmlyaWRpcyg2KSksIDAuNykpXG5cbmdncGxvdChwcm9wX2hpZ2hfdmlydXNlc19sb25nW3Byb3BfaGlnaF92aXJ1c2VzX2xvbmckdnM9PTEsXSwgYWVzKHg9YXMuZmFjdG9yKG51bXJ1bGVzKSwgeT1wcm9wb3J0aW9uLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb2xvcj1hcy5mYWN0b3IobnVtcnVsZXMpLCBmaWxsPWFzLmZhY3RvcihudW1ydWxlcykpKSArXG4gIGdlb21fYm94cGxvdCgpICtcbiAgI2dlb21fcG9pbnQoKSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCwgYW5nbGUgPSA5MCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTYpLFxuICApICtcbiAgeGxhYihcIk51bWJlciBvZiBSdWxlc1wiKSArXG4gIHlsYWIoXCJQcm9wb3J0aW9uIG9mIENvbnRpZ3MgaW4gQ29tbW9uXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAwLjMpKSArXG4gIHNjYWxlX2NvbG9yX21hbnVhbChuYW1lPVwiXCIsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShyZXYodmlyaWRpcyg2KSksIDAuNykpXG5cbmdncGxvdChwcm9wX2hpZ2hfdmlydXNlc19sb25nW3Byb3BfaGlnaF92aXJ1c2VzX2xvbmckdmI9PTEsXSwgYWVzKHg9YXMuZmFjdG9yKG51bXJ1bGVzKSwgeT1wcm9wb3J0aW9uLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBjb2xvcj1hcy5mYWN0b3IobnVtcnVsZXMpLCBmaWxsPWFzLmZhY3RvcihudW1ydWxlcykpKSArXG4gIGdlb21fYm94cGxvdCgpICtcbiAgI2dlb21fcG9pbnQoKSArXG4gIHRoZW1lX2xpZ2h0KCkgK1xuICB0aGVtZShcbiAgICBwYW5lbC5ncmlkLm1ham9yLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgcGFuZWwuYm9yZGVyID0gZWxlbWVudF9ibGFuaygpLFxuICAgIGF4aXMudGlja3MueSA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBsZWdlbmQucG9zaXRpb24gPSBcImJvdHRvbVwiLFxuICAgIGF4aXMudGV4dC55PWVsZW1lbnRfdGV4dChzaXplPTE0KSxcbiAgICBheGlzLnRleHQueD1lbGVtZW50X3RleHQoc2l6ZT0xNCwgYW5nbGUgPSA5MCksXG4gICAgbGVnZW5kLnRleHQ9ZWxlbWVudF90ZXh0KHNpemU9MTIpLFxuICAgIGF4aXMudGl0bGU9ZWxlbWVudF90ZXh0KHNpemU9MTYpLFxuICApICtcbiAgeGxhYihcIk51bWJlciBvZiBSdWxlc1wiKSArXG4gIHlsYWIoXCJQcm9wb3J0aW9uIG9mIENvbnRpZ3MgaW4gQ29tbW9uXCIpICtcbiAgc2NhbGVfZmlsbF9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAwLjMpKSArXG4gIHNjYWxlX2NvbG9yX21hbnVhbChuYW1lPVwiXCIsXG4gICAgICAgICAgICAgICAgICAgICB2YWx1ZXMgPSBhbHBoYShyZXYodmlyaWRpcyg2KSksIDAuNykpXG5cbmdncGxvdChwcm9wX2hpZ2hfdmlydXNlc19sb25nW3Byb3BfaGlnaF92aXJ1c2VzX2xvbmckZHZmPT0xLF0sIGFlcyh4PWFzLmZhY3RvcihudW1ydWxlcyksIHk9cHJvcG9ydGlvbiwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29sb3I9YXMuZmFjdG9yKG51bXJ1bGVzKSwgZmlsbD1hcy5mYWN0b3IobnVtcnVsZXMpKSkgK1xuICBnZW9tX2JveHBsb3QoKSArXG4gICNnZW9tX3BvaW50KCkgK1xuICB0aGVtZV9saWdodCgpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuZ3JpZC5tYWpvci55ID0gZWxlbWVudF9ibGFuaygpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBheGlzLnRpY2tzLnkgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgbGVnZW5kLnBvc2l0aW9uID0gXCJib3R0b21cIixcbiAgICBheGlzLnRleHQueT1lbGVtZW50X3RleHQoc2l6ZT0xNCksXG4gICAgYXhpcy50ZXh0Lng9ZWxlbWVudF90ZXh0KHNpemU9MTQsIGFuZ2xlID0gOTApLFxuICAgIGxlZ2VuZC50ZXh0PWVsZW1lbnRfdGV4dChzaXplPTEyKSxcbiAgICBheGlzLnRpdGxlPWVsZW1lbnRfdGV4dChzaXplPTE2KSxcbiAgKSArXG4gIHhsYWIoXCJOdW1iZXIgb2YgUnVsZXNcIikgK1xuICB5bGFiKFwiUHJvcG9ydGlvbiBvZiBDb250aWdzIGluIENvbW1vblwiKSArXG4gIHNjYWxlX2ZpbGxfbWFudWFsKG5hbWU9XCJcIixcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IGFscGhhKHJldih2aXJpZGlzKDYpKSwgMC4zKSkgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwobmFtZT1cIlwiLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gYWxwaGEocmV2KHZpcmlkaXMoNikpLCAwLjcpKVxuYGBgIn0= -->

```r
ggplot(prop_high_viruses_long[prop_high_viruses_long$vs2==1,], aes(x=as.factor(numrules), y=proportion, 
                                  color=as.factor(numrules), fill=as.factor(numrules))) +
  geom_boxplot() +
  #geom_point() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Number of Rules") +
  ylab("Proportion of Contigs in Common") +
  scale_fill_manual(name="",
                     values = alpha(rev(viridis(6)), 0.3)) +
  scale_color_manual(name="",
                     values = alpha(rev(viridis(6)), 0.7))

ggplot(prop_high_viruses_long[prop_high_viruses_long$vs==1,], aes(x=as.factor(numrules), y=proportion, 
                                  color=as.factor(numrules), fill=as.factor(numrules))) +
  geom_boxplot() +
  #geom_point() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Number of Rules") +
  ylab("Proportion of Contigs in Common") +
  scale_fill_manual(name="",
                     values = alpha(rev(viridis(6)), 0.3)) +
  scale_color_manual(name="",
                     values = alpha(rev(viridis(6)), 0.7))

ggplot(prop_high_viruses_long[prop_high_viruses_long$vb==1,], aes(x=as.factor(numrules), y=proportion, 
                                  color=as.factor(numrules), fill=as.factor(numrules))) +
  geom_boxplot() +
  #geom_point() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Number of Rules") +
  ylab("Proportion of Contigs in Common") +
  scale_fill_manual(name="",
                     values = alpha(rev(viridis(6)), 0.3)) +
  scale_color_manual(name="",
                     values = alpha(rev(viridis(6)), 0.7))

ggplot(prop_high_viruses_long[prop_high_viruses_long$dvf==1,], aes(x=as.factor(numrules), y=proportion, 
                                  color=as.factor(numrules), fill=as.factor(numrules))) +
  geom_boxplot() +
  #geom_point() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 90),
    legend.text=element_text(size=12),
    axis.title=element_text(size=16),
  ) +
  xlab("Number of Rules") +
  ylab("Proportion of Contigs in Common") +
  scale_fill_manual(name="",
                     values = alpha(rev(viridis(6)), 0.3)) +
  scale_color_manual(name="",
                     values = alpha(rev(viridis(6)), 0.7))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


# Comparing RefSeq vs VS2 database
VS2 sequences used: (/scratch/kwigg_root/kwigg/hegartyb/VSTE/non-refseq-genomes/NCLDV/env/NCLDV_env.fna /scratch/kwigg_root/kwigg/hegartyb/VSTE/non-refseq-genomes/lavidaviridae/env/HQ_virophages.fna /scratch/kwigg_root/kwigg/hegartyb/VSTE/non-refseq-genomes/dsDNAphage/env/* /scratch/kwigg_root/kwigg/hegartyb/VSTE/non-refseq-genomes/dsDNAphage/provirus/VirSorterCuratedDataset_HQ_contigs.fna)


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxuYW5uX2NvbCRudW1ydWxlcyA8LSByb3dTdW1zKGFubl9jb2xbLDE6Nl0pXG5hbm5fY29sJG51bXRvb2xzIDwtIGFubl9jb2wkbnVtcnVsZXNcblxuZm9yIChpIGluIDE6bnJvdyhhbm5fY29sKSkge1xuICBpZiAoYW5uX2NvbCR0dltpXT09MSkge1xuICAgIGFubl9jb2wkbnVtdG9vbHNbaV0gPC0gYW5uX2NvbCRudW10b29sc1tpXSArIDNcbiAgfVxuICBcbiAgaWYgKGFubl9jb2wkdG52W2ldPT0xKSB7XG4gICAgYW5uX2NvbCRudW10b29sc1tpXSA8LSBhbm5fY29sJG51bXRvb2xzW2ldICsgNVxuICB9XG59XG5cbmFubl9jb2wkbnVtdG9vbHNbYW5uX2NvbCRudW10b29scz42XSA8LSA2XG5cbmFubl9jb2wkdG9vbCA8LSByb3duYW1lcyhhbm5fY29sKVxuYGBgXG5gYGAifQ== -->

```r
```r
ann_col$numrules <- rowSums(ann_col[,1:6])
ann_col$numtools <- ann_col$numrules

for (i in 1:nrow(ann_col)) {
  if (ann_col$tv[i]==1) {
    ann_col$numtools[i] <- ann_col$numtools[i] + 3
  }
  
  if (ann_col$tnv[i]==1) {
    ann_col$numtools[i] <- ann_col$numtools[i] + 5
  }
}

ann_col$numtools[ann_col$numtools>6] <- 6

ann_col$tool <- rownames(ann_col)
```
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGFibGUocmVmc2VxX3ZpciRpZCAlaW4lIHZpcnVzZXMkY29udGlnKVxuYGBgIn0= -->

```r
table(refseq_vir$id %in% viruses$contig)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlc19yZWZzZXFfdmlyIDwtIHZpcnVzZXNbdmlydXNlcyRjb250aWcgJWluJSByZWZzZXFfdmlyJGlkLF1cbnZpcnVzZXNfbm90X3JlZnNlcV92aXIgPC0gdmlydXNlc1shKHZpcnVzZXMkY29udGlnICVpbiUgcmVmc2VxX3ZpciRpZCkgJiB2aXJ1c2VzJHNlcXR5cGU9PVwidmlydXNcIixdXG5gYGAifQ== -->

```r
viruses_refseq_vir <- viruses[viruses$contig %in% refseq_vir$id,]
viruses_not_refseq_vir <- viruses[!(viruses$contig %in% refseq_vir$id) & viruses$seqtype=="virus",]
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxucHJvcF9oaWdoX3ZpcnVzZXNfbG9uZyA8LSBhcy5kYXRhLmZyYW1lKHByb3BfaGlnaF92aXJ1c2VzKVxucHJvcF9oaWdoX3ZpcnVzZXNfbG9uZyR0b29sIDwtIGNvbWJvc19saXN0JHRvb2xjb21ib1xuXG5wcm9wX2hpZ2hfdmlydXNlc19sb25nIDwtIGFzLmRhdGEuZnJhbWUocHJvcF9oaWdoX3ZpcnVzZXNfbG9uZykgJT4lIFxuICBwaXZvdF9sb25nZXIoIXRvb2wsXG4gICAgICAgICAgICAgICB2YWx1ZXNfdG89XFxwcm9wb3J0aW9uXFwsXG4gICAgICAgICAgICAgICBuYW1lc190bz1cXHRvb2wyXFwpICU+JVxuICBpbm5lcl9qb2luKGFubl9jb2wsIGJ5PWMoXFx0b29sXFwpKVxuXG5cbmBgYFxuYGBgIn0= -->

```r
```r
prop_high_viruses_long <- as.data.frame(prop_high_viruses)
prop_high_viruses_long$tool <- combos_list$toolcombo

prop_high_viruses_long <- as.data.frame(prop_high_viruses_long) %>% 
  pivot_longer(!tool,
               values_to=\proportion\,
               names_to=\tool2\) %>%
  inner_join(ann_col, by=c(\tool\))

```
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudmlydXNlc19yZWZzZXFfdmlyJGtlZXBfc2NvcmVfaGlnaF9NQ0MgPC0gZ2V0dGluZ192aXJhbF9zZXRfMSh2aXJ1c2VzX3JlZnNlcV92aXIsIGluY2x1ZGVfZGVlcHZpcmZpbmRlciA9IEYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aWJyYW50ID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpcnNvcnRlcjIgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX3ZpcmFsID0gVCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3R1bmluZ19ub3RfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyID0gRilcblxudmlydXNlc19ub3RfcmVmc2VxX3ZpciRrZWVwX3Njb3JlX2hpZ2hfTUNDIDwtIGdldHRpbmdfdmlyYWxfc2V0XzEodmlydXNlc19ub3RfcmVmc2VxX3ZpciwgaW5jbHVkZV9kZWVwdmlyZmluZGVyID0gRixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbmNsdWRlX3ZpYnJhbnQgPSBGLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdmlyc29ydGVyMiA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV90dW5pbmdfdmlyYWwgPSBULFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGluY2x1ZGVfdHVuaW5nX25vdF92aXJhbCA9IFQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5jbHVkZV92aXJzb3J0ZXIgPSBGKVxuYGBgIn0= -->

```r
viruses_refseq_vir$keep_score_high_MCC <- getting_viral_set_1(viruses_refseq_vir, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)

viruses_not_refseq_vir$keep_score_high_MCC <- getting_viral_set_1(viruses_not_refseq_vir, include_deepvirfinder = F,
                                              include_vibrant = F,
                                              include_virsorter2 = T,
                                              include_tuning_viral = T,
                                              include_tuning_not_viral = T,
                                              include_virsorter = F)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuc3VtKHZpcnVzZXNfcmVmc2VxX3ZpciRrZWVwX3Njb3JlX2hpZ2hfTUNDPj0xKS9ucm93KHZpcnVzZXNfcmVmc2VxX3ZpcilcbmBgYCJ9 -->

```r
sum(viruses_refseq_vir$keep_score_high_MCC>=1)/nrow(viruses_refseq_vir)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuc3VtKHZpcnVzZXNfbm90X3JlZnNlcV92aXIka2VlcF9zY29yZV9oaWdoX01DQz49MSkvbnJvdyh2aXJ1c2VzX25vdF9yZWZzZXFfdmlyKVxuYGBgIn0= -->

```r
sum(viruses_not_refseq_vir$keep_score_high_MCC>=1)/nrow(viruses_not_refseq_vir)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->






<!-- rnb-text-end -->

