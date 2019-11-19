---
title: 김정주의 블로그
filename: index.md
---
## 스터디

[로짓이란?](study/2019-11-19-logit.md)

          {% for page in site.pages %}
          <p class="view">
            <a href="{{ page.filename }}">{{ page.title }}</a>
          </p>
          {% endfor %}
