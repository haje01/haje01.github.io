---
layout: template
filename: 2019-11-19-github-page.md
---

# Jekyll 없이 개인용 Github 페이지 만들기

Jekyll 설치 없이 내 Github 페이지에 개인용 블로그를 만든 기록.

1. 먼저 `GITHUB_ID.github.io` 식으로 저장소를 만듦
2. 다음 페이지에서 `Settings` 클릭
3. 다음 페이지의 `GitHub Pages` 섹션에서 `Choose a theme`을 클릭
4. 원하는 테마 선택 후 `Select theme` 클릭
5. 기본 `index.md` 파일의 내용을 아래와 같이 채우고 커밋

```
---
title: 나의 노트
layout: template
filename: index.md
---
```
6. 저장소를 git으로 로컬에 클론
7. 클론한 폴더 아래 `_layout` 폴더를 만들고, 아래와 같이 `template.html` 파일을 만듦.

```html
{% raw %}
<!DOCTYPE html>
<html lang="en-US">
  <head>
    <meta charset="UTF-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

<!-- Begin Jekyll SEO tag v2.5.0 -->
<title>나의 노트 | GITHUB-ID.github.io</title>
<meta name="generator" content="Jekyll v3.8.5" />
<meta property="og:title" content="나의 노트" />
<meta property="og:locale" content="en_US" />
<link rel="canonical" href="https://GITHUB-ID.github.io/" />
<meta property="og:url" content="https://GITHUB-ID.github.io/" />
<meta property="og:site_name" content="GITHUB-ID.github.io" />
<script type="application/ld+json">
{"@type":"WebSite","url":"https://GITHUB-ID.github.io/","name":"나의 노트","headline":"나의 노트","@context":"http://schema.org"}</script>
<!-- End Jekyll SEO tag -->

    <link rel="stylesheet" href="/assets/css/style.css?v=01e6290648d6409b0c7f076e8788b0cbc74c3e34">

    <!-- MathJax 설정 -->
    <script type="text/x-mathjax-config">
      MathJax.Hub.Config({
        extensions: ["tex2jax.js"],
        jax: ["input/TeX", "output/HTML-CSS"],
        tex2jax: {
          inlineMath: [ ['$','$'], ["\\(","\\)"] ],
          displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
          processEscapes: true
        },
        "HTML-CSS": { availableFonts: ["TeX"] }
      });
    </script>
    <script type="text/javascript" async
    src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML">
    </script>
    <!--[if lt IE 9]>
    <script src="//cdnjs.cloudflare.com/ajax/libs/html5shiv/3.7.3/html5shiv.min.js"></script>
    <![endif]-->
  </head>
  <body>
    <div class="wrapper">
      <header>
        <h1><a href="https://GITHUB-ID.github.io/">나의 노트</a></h1>

        <!-- 글의 분류 -->
        <h2>스터디</h2>
          {% for page in site.pages %}
            {% if page.dir == "/study/" %}
            <p class="view">
              <a href="https://GITHUB-ID.github.io/study/{{ page.filename | remove: ".md" }}.html">{{ page.title }}</a>
            </p>
            {% endif %}
          {% endfor %}

      </header>
      <section>
{{ content }}
      </section>
      <footer>

        <p><small>Hosted on GitHub Pages &mdash; Theme by <a href="https://github.com/orderedlist">orderedlist</a></small></p>
      </footer>
    </div>
    <script src="/assets/js/scale.fix.js"></script>

  </body>
</html>
{% endraw %}
```
8. 이제 필요에 따라 글을 분류할 하위 폴더들을 만들어줌. 예를 들어 공부 노트를 위한 `study` 폴더.
9. 분류에 맞는 폴더에 Markdown 형식으로 글을 작성하고, `study/2019-11-19-SUBJECT-NAME.md` 식으로 저장
10. 단, 아래와 같은 Jekyll용 메타 정보가 제일 위에 와야함.
```
---
layout: template
filename: 2019-11-19-SUBJECT-NAME
---
```
10. 글을 커밋하고, 웹브라우저에서 `GITHUB-ID.github.io`를 방문해 잘 나오는지 확인.

## 참고 링크
* https://phuston.github.io/patrickandfrantonarethebestninjas/howto
