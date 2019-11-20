---
layout: template
title: Jekyll 없이 개인용 Github 페이지 만들기
tags: git
---

Jekyll 설치 없이 내 Github 페이지에 개인용 블로그를 만든 기록.

## 기본 설정

* 먼저 `GITHUB_ID.github.io` 식으로 저장소를 만듦
* 다음 페이지에서 `Settings` 클릭
* 다음 페이지의 `GitHub Pages` 섹션에서 `Choose a theme`을 클릭
* 원하는 테마 선택 후(여기에서는 minimal 테마 기준으로 설명) `Select theme` 클릭
* 기본 `index.md` 파일의 내용을 아래와 같이 채우고 커밋

```
---
title: 나의 노트
layout: template
filename: index.md
---
```
* 저장소를 git으로 로컬에 클론
* 기본 폴더(클론한 폴더) 아래 `_layout` 폴더를 만들고, 아래와 같이 `template.html` 파일을 만듦.

```html
{% raw %}
<!DOCTYPE html>
<html lang="en-US">
  <head>
    <meta charset="UTF-8">`
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

<!-- Begin Jekyll SEO tag v2.5.0 -->
<title>김정주의 노트 | haje01.github.io</title>
<meta name="generator" content="Jekyll v3.8.5" />
<meta property="og:title" content="김정주의 노트" />
<meta property="og:locale" content="en_US" />
<link rel="canonical" href="https://haje01.github.io/" />
<meta property="og:url" content="https://haje01.github.io/" />
<meta property="og:site_name" content="haje01.github.io" />
<script type="application/ld+json">
{"@type":"WebSite","url":"https://haje01.github.io/","name":"김정주의 노트","headline":"김정주의 노트","@context":"http://schema.org"}</script>
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

    <!-- 태그 수집 -->
    {% assign rawtags = "" %}
    {% for post in site.posts %}
      {% assign ttags = post.tags | join:'|' | append:'|' %}
      {% assign rawtags = rawtags | append:ttags %}
    {% endfor %}
    {% assign rawtags = rawtags | split:'|' | sort %}

    {% assign all_tags = "" %}
    {% for tag in rawtags %}
      {% if tag != "" %}
        {% if all_tags == "" %}
          {% assign all_tags = tag | split:'|' %}
        {% endif %}
        {% unless all_tags contains tag %}
          {% assign all_tags = all_tags | join:'|' | append:'|' | append:tag | split:'|' %}
        {% endunless %}
      {% endif %}
    {% endfor %}

    <!--[if lt IE 9]>
    <script src="//cdnjs.cloudflare.com/ajax/libs/html5shiv/3.7.3/html5shiv.min.js"></script>
    <![endif]-->
  </head>
  <body>
    <div class="wrapper">
      <header>
        <h1><a href="https://haje01.github.io/">김정주의 노트</a></h1>

          {% for post in site.posts %}
            <p class="view">
              <a href="{{ post.url }}">{{ post.title }}</a>
            </p>
          {% endfor %}

      </header>
      <section>
      {% if page.tags and page.title | split:' ' | first != "Tag:" %}
      <span style="color: purple">[
          {% for tag in page.tags %}
            {% capture tag_name %}{{ tag }}{% endcapture %}
            <a href="/tag/{{ tag_name }}"><code class="highlighter-rouge"><nobr>{{ tag_name }}</nobr></code>&nbsp;</a>
          {% endfor %}
        ]</span>
      {% endif %}


      {% if page.title %}
      <h1>{{ page.title }}</h1>
      {% endif %}

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

## 태그 이용하기
포스트별로 하나 이상의 태그를 붙여 관리할 수 있다. 먼저 몇 가지 작업이 필요하다.

* `_layouts` 폴더 아래에 `tagpage.html` 파일을 아래와 같은 내용으로 작성한다:

```html
{% raw %}
---
layout: template
---
<div class="post">
<ul>
{% for post in site.tags[page.tag] %}
  <li><a href="{{ post.url }}">{{ post.title }}</a> ({{ post.date | date: "%Y-%m-%d" }})<br>
    {{ post.description }}
  </li>
{% endfor %}
</ul>
</div>
<hr>
{% endraw %}
```

* `tag` 폴더를 만들고, 사용할 태그 타입별로 아래와 같은 파일을 `TAG-NAME.md` 식으로 저장하여야 한다.

```
---
layout: tagpage
title: "Tag: TAG-NAME"
tag: TAG-NAME
---
```

## 글의 작성

* 이제 기본 폴더 아래 `_posts` 폴더를 만들고, 거기에 Markdown 형식으로 글을 작성한다. (예: `_posts/2019-11-19-SIMPLE-SUBJECT.md`)
* 단, 아래와 같은 Jekyll용 메타 정보가 제일 위에 와야한다.

```
---
layout: template
title: POST-TITLE
tag: TAG-NAME
---
```

* 글을 커밋하고, 웹브라우저에서 `GITHUB-ID.github.io`를 방문해 잘 나오는지 확인.

정상적으로 동작하면, 각 글의 제목 위에 태그를 볼 수 있으며, 태그명을 누르면 그 태그에 속한 모든 글의 리스트를 볼 수 있다.

## 참고 링크
* https://phuston.github.io/patrickandfrantonarethebestninjas/howto
* https://longqian.me/2017/02/09/github-jekyll-tag/