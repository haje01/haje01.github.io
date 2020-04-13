---
layout: post
title: 스네이크메이크 (Snakemake) 튜토리얼
description:
date: 2020-04-02
tags: [knowhow, draft]
---

이 글에서는 파이썬용 워크플로우 관리 시스템인 Snakemake (스네이크메이크) 에 대해 살펴보겠다. Snakemake 버전 5.14 와 유닉스 계열 (Linux/macOS) OS 관점에서 설명하지만, 경로 구분자 등 몇 가지만 고려하면 윈도우에서도 적용에 무리가 없을 것이다.

## 워크플로우 관리가 필요한 이유

요즘은 데이터 과학과 인공지능이 각광받는 시대이지만, 실상 그것을 위한 준비 과정, 즉 데이터 엔지니어링 과정의 어려움은 간과되기 쉽다. 간단한 데이터 처리라면 그다지 복잡할 것이 없겠지만, 현업에서 높은 결과를 내기 위해 다양한 소스의 데이터를 처리(ETL) 하고, 그것에서 고도화된 피처 엔지니어링을 수행하다 보면 어려움이 많다:

* 다양한 데이터 소스의 의존성 파악이 힘듦 - 작업에 필요한 데이터가 적절히 준비되어 있나?
* 반복적인 데이터 처리작업으로 인한 비효율성 - 필요하지 않은데 했던 처리를 반복하는 경우가 있음
* 처리 코드와 의존 관계가 섞여 재활용성이 낮아지고 가독성이 떨어짐 - 데이터 처리와 의존성을 해결을 분리할 수 없을까?

![데이터 의존 관계과 처리를 함께 기술](/assets/2020-04-07-11-56-29.png)

이와 같은 문제를 해결하기 위해 데이터 처리용 다양한 워크플로우 관리 시스템 (Workflow Management System, WMS) 이 활용되고 있다. 대부분의 워크플로우 매니저는 유향 비순환 그래프 (Directed Acyclic Graph, DAG) 방식으로 구현이 되어 있다. 이 그래프에서 각 노드는 작업을 의미하고, 자기 작업을 수행하기 위해 필요한 의존 작업의 노드에 연결되어 있다. 타겟 결과물이 요청되면 그것을 생성하기 위해 필요한 의존 작업을 재귀적으로 거슬러 올라가며 호출하는 방식이다.

[AirFlow](https://airflow.apache.org), [Luigi](https://github.com/spotify/luigi) 등 다양한 WMS 가 쓰이고 이는데, 이 글에서는 비교적 많이 알려지지 않은 Snakemake (스네이크메이크) 에 대해서 살펴보겠다. 내가 생각하는 Snakemake 의 장점은:

* GNU Make 형식의 컴팩트한 구문
* 데이터 의존관계 파악이 쉽다.
* 처리 코드와 의존 관계 명세를 나누어 재활용성이 높다.
* Python으로 구현되어 호환성이 좋다.

위와 같은 장점에도, Snakemake 가 널리 쓰이지 않는 이유는, Snakemake의 개발자가 Bio Informatics 쪽을 전문으로 하고 있어서, 관련 예제가 전부 바이오 쪽의 툴과 데이터를 이용하기 때문인 것 같다. 이에 이 예제에서는 가급적 일반적인 예로 설명하려 한다.

![데이터 의존 관계와 처리를 별도로 기술](/assets/2020-04-07-12-13-53.png)

## 빠르게 시작하기

요즘처럼 통합개발환경 (IDE) 가 없었던 오래전에는, 프로그램을 빌드하기 위해서 여러가지 소스 코드 파일과 중간 산출물 파일들을 직접 관리해 주어야 했는데 이것은 굉장히 복잡한 과정이었다. 이를 해결하기 위해서 사용된 것이 GNU Make인데, Snakemake는 그것의 철학을 이어받은 WMS 라고 할 수 있겠다. 간단한 데이터 처리 예롤 살펴보며 Snakemake의 기본 개념을 확인해 보자.

현재 디렉토리가 다음과 같이 구성되어 있다고 하자:

```
    data/
        A.txt
        B.txt
        C.txt
```

파일의 내용은 다음과 같다.

`A.txt`
```
Hello.
```

`B.txt`
```
This is your captain speaking.
```

`C.txt`
```
Welcome aboard.
```

지금부터 하려고 하는 것은, 각 파일에 들어있는 단어의 수를 세고, 그것을 그래프로 만드는 일이다. 이를 위해서는 `Snakefile` 이라는 이름으로 워크플로우 파일이 필요하다. Snakemake 는 이 워크플로우 파일의 명세대로 동작하게 된다.

일반적으로 워크플로우는 하나 이상의 작업 노드들이 연결된 DAG 형식으로 되어 있는데, Snakemake 에서는 이 노드를 규칙(Rule) 이라고 부르며, `rule` 이라는 키워드로 선언한다. 규칙은 기본적으로 다음과 같은 구조를 가진다:

```
rule 규칙 이름:
    입력 파트:
    출력 파트:
    실행 파트:
```

파트별 구체적인 키워드는 다음과 같다.

```
rule RULE_NAME:
    input:
    output:
    run, script, shell 중 택일:
```

* `rule` 다음에는 규칙의 이름이 온다.
* `input:` 다음에는 이 규칙이 사용하는 입력 파일 이름이나 패턴이 온다. 하나 이상인 경우 `,` 를 사용하여 구분해야 한다.
* `output` 다음에는 이 규칙의 실행 결과로 출력되는 파일 이름이나 패턴이 온다. 하나 이상인 경우 `,` 를 사용하여 구분해야 한다.
* 실행 파트에는 다음과 같은 키워드가 온다
  * `shell` - 실행할 쉘스크립트를 기술
  * `run` - 파이썬 스크립트를 기술
  * `script` - 외부 파이썬 스크립트 파일명을 지정

> `Snakefile` 은 기본적으로 파이썬 파일과 같다. 즉, 파이썬 모듈을 임포트하거나 파이썬 명령어를 사용할 수 있다. 단, 위에서 처럼 **규칙 (Rule)** 구획은 다른 신택스를 사용한다.

### 파일내 단어수 세기

이제 다음과 같은 내용으로 `Snakefile` 라는 파일을 현재 디렉토리에 만들어 보자:

```python
rule count:
    """파일내 단어수 세기."""
    input:
        "data/{filename}.txt
    output:
        "temp/wc_{filename}.txt
    shell:
        "wc -w {input}
```

> 입력 및 출력에서 `{filename}` 부분은 **와일드카드** 로 불리는데, 패턴이 매칭되는 모든 파일에서 대체되는 변수이다.

> `wc -w` 는 파일에서 단어를 세는 Unix 계열 쉘 명령어다. 예를 들어 `data/A.txt` 파일의 단어를 세려면 다음과 같이 한다:
> ```
> $ wc -w data/A.txt
>        1 data/A.txt
> ```

위 예는 `count` 라는 규칙이 `data/` 디렉토리 아래의 `[파일명].txt` 에 매치되는 모든 파일을 입력으로 하여 `shell` 에 명기된 명령을 수행한 후, 그 결과를 `temp/wc_[파일명].txt` 파일로 저장한다는 뜻이다.

예로, 이 Snakefile이 있는 디렉토리에서 다음과 같이 실행해보자:

```
$ snakemake --cores 1 temp/wc_A.txt
```

위 명령은 `temp/wc_A.txt` 를 생성하라는 것인데, 이를 위해서는 출력이 `temp/wc_A.txt` 에 매칭되는 규칙이 필요한데, 그 규칙은 `count` 이다. 그것에 따라 `A.txt` 파일의 단어수를 `wc -w` 명령으로 세어 `result` 디렉토리의 `wc_A.txt` 파일에 저장한다. 다음과 같이 결과를 확인해보자:

```
$ cat temp/wc_A.txt

       1 data/A.txt
```

이것은 `A.txt` 파일에는 하나의 단어가 있다는 뜻이다. 같은 식으로 `B.txt` 와 `C.txt` 에 대해서도 단어 수를 구할 수 있다.

### 개발 단어수 파일 병합

다음 단계는 각 파일별 단어 수를 하나의 파일로 결합하고 `.csv` 형식으로 저장하는 것이다. `Snakefile` 에 아래처럼 `concat` 규칙을 추가한다:

```python
rule concat:
    """개별 단어수 파일을 병합."""
    input:
        expand("temp/wc_{filename}.txt", filename=['A', 'B', 'C'])
    output:
        "temp/wc_all.csv"
    script:
        "concat.py"
```

위 규칙에서는 모든 단어 수 파일을 읽어 하나의 파일로 출력한다. 즉, `count` 규칙은 입력과 출력 파일이 일대일의 관계였지만, `concat` 규칙은 다대일의 관계인 것이다. 이를 위해서 `input` 에서 모든 단어 수 파일을 명시해주어야 하는데, 이를 위해 `expand` 함수를 사용하고 있다.

`expand` 는 매개변수 리스트의 각 항목으로 특정 패턴을 확장한다. 위의 경우, `input` 의 값은 아래와 같다:

```
"temp/wc_A.txt", "temp/wc_B.txt", "temp/wc_C.txt",
```

> `expand` 의 인자 내 `{filename}` 은 와일드카드가 아니고, 매개변수 리스트를 통해 확장될 플레이스 홀더이다. `expand` 인자에서 와일드카드를 쓰려면 이중 중괄호 `{{}}` 를 사용해야 한다.

`Snakefile` 은 기본적으로 파이썬 파일이기에, 다음과 같이 `FILENAMES` 변수를 이용할 수도 있다.
```python
FILENAMES = ['A', 'B', 'C']

rule concat:
    """개별 단어수 파일을 병합."""
    input:
        expand("temp/wc_{filename}.txt", filename=FILENAMES)
    output:
        "temp/wc_all.csv"
    script:
        "concat.py"
```


또 이 규칙에서는 세 파일 결합하기 위해 `concat.py` 이라는 파이썬 스크립트 파일을 실행하여 `temp/wc_all.csv` 를 생성하는데, 그것은 내용은 다음과 같다:

```python
import re
import pandas as pd

# wc -w 결과 파싱용 정규식
PTRN = re.compile(r'\s+(\d+)\sdata/(.*)$')

# 출력용 csv 파일
with open(snakemake.output[0], 'wt') as f:
    f.write('fname, count\n')
    # snakefile에 명시된 모든 입력 파일에 대해서
    for fn in snakemake.input:
        # 파싱하고 출력
        line = open(fn, 'rt').read()
        cnt, fn = PTRN.search(line).groups()
        f.write(f'{fn}, {cnt}\n')
```

Snakemake 를 통해서 실행되는 파이썬 스크립트에는 `snakemake` 라는 객체가 기본적으로 생성되며, 여기에 Snakemake 관련 다양한 정보를 속성으로 갖는다. 위 스크립트에서는 `snakemake` 의 `concat` 규칙에서 명시된 입력 파일을 얻어 오기 위해서 `snakemake.input` 을, 출력 파일을 얻어 오기 위해서 `snakemake.output` 을 이용하고 있다. `Snakefile` 의 `input` 이나 `output` 에는 하나 이상의 파일이 올 수 있기에, 배열처럼 인덱스를 사용하고 있다.

> 스크립트에서 파일 입출력을 작성할 때, 특정 파일명이 아닌 `snakemake.input` 이나 `snakemake.output` 을 사용해야 한다는 점에 다시 한 번 주의하자.

이제 다음처럼 실행하면:

```
$ snakemake --cores 1 temp/wc_all.csv
```

`A.txt`, `B.txt`, `C.txt` 세 파일의 단어 수가 아래처럼 하나의 `.csv` 파일에 출력된다.

```
$ cat temp/wc_all.csv

fname, count
A.txt, 1
B.txt, 5
C.txt, 2
```

### 그래프 그리기

위에서 생성된 `.csv` 를 읽어, 그래프를 그려보자. `snakefile` 에 아래와 같은 규칙을 추가한다:

```python
rule plot:
    input:
        "temp/wc_all.csv"
    output:
        "temp/wc_all.png"
    script:
        "plot.py"
```

`wc_all.csv` 파일을 읽어 그래프를 그린 결과를 `wc_all.png` 에 저장한다. 이를 위한 스크립트파일 `plot.py` 의 내용은 아래와 같다:

```python
import pandas as pd

df = pd.read_csv(snakemake.input[0], index_col='fname')
plot = df.plot(kind='bar', rot=45)
fig = plot.get_figure()
fig.tight_layout()
fig.savefig(snakemake.output[0])
```

이제 실행하면 멋진(?) 그래프가 만들어진다.

```
$ snakemake --cores 1 temp/wc_all.png
```

![단어 그래프](/assets/2020-04-08-14-27-30.png)

참고로 이제 `snakefile` 의 내용은 아래와 같을 것이다.

```python
FILENAMES = ['A', 'B', 'C']


rule all:
    input:
        "temp/wc_all.png"

rule count:
    input:
        "data/{filename}.txt"
    output:
        "temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    input:
        expand('temp/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/wc_all.csv"
    script:
        "concat.py"

rule plot:
    input:
        "temp/wc_all.csv"
    output:
        "temp/wc_all.png"
    script:
        "plot.py"
```

### 의존 관계 톮아보기

지금까지 각 단계별로 규칙을 만들고, 입력 데이터나 한 규칙의 결과물을 다음 규칙이 이용해 결과물을 내는 것을 살펴보았다. 우리는 단계별로 진행을 했기에 중간 결과물들이 차례되로 생성되어 있었지만, 중간 결과물이 없을 때 바로 타겟 파일 생성을 요구하면 DAG 에 따라 아래와 같이 진행된다.

1. 명시적 타겟 파일 없이 Snakemake 실행
2. `all` 규칙에서 `temp/wc_all.png` 파일 생성을 요구
3. `plot` 규칙에서 `temp/wc_all.csv` 파일 생성을 요구
4. `concat` 규칙에서 `temp/wc_A.txt`, `temp/wc_B.txt`, `temp/wc_C.txt` 파일 생성을 요구
5. `count` 규칙에서 `data/A.txt`, `data/B.txt`, `data/C.txt` 의 단어를 세고 결과를 각각 `temp/wc_A.txt`, `temp/wc_B.txt`, `temp/wc_C.txt` 에 저장
6. `concat` 규칙에서 `temp/wc_all.csv` 파일을 생성
7. `plot` 규칙에서 `temp/wc_all.png` 파일을 생성

이러한 의존 관계를 통한 워크플로우 관리는 미리 순차적으로 작업을 하는 것이 아니라, 요청된 타겟 파일의 생성에 필요한 의존성을 거꾸로 쫓아 재귀적으로 진행된다. 일종의 **On-Demand** 또는 **Lazy Evaluation** 으로 볼 수 있겠다. 이 방식은 필요할 때 필요한 작업만 진행되기에 준비 작업이 적고, 앞 규칙의 저정된 결과를 그대로 이용하기에, 부분 규칙을 여러번 재작업해도 최단 시간에 결과를 확인할 수 있다.

예를 들어 앞의 `plot` 규칙에서, 좀 더 미려한 그래프를 그리기 위해서 코드를 여러번 수정해야 한다면, (Jupyter 노트북처럼 인터랙티브한 환경이 아닌 경우) 매번 실행하여 결과를 확인해야한다. 이 과정에서 앞에서 실행했던 `count` 나 `concat` 이 불필요하게 실행될 수 있겠다. 물론, 각 단계에서 중간 결과를 저장하고, 중간 결과가 있는 경우 그 것을 이용하도록 코드를 수정할 수도 있겠으나, 아무래도 번거롭고 코드가 지저분해지게 된다. 예에서 처럼 규칙을 통해 의존 관계를 선언하면 앞 단계 규칙의 결과, 즉 `temp/wc_all.csv`가 있는 경우 빠르게 그래프 그리는 코드만 실행하면 되기에 효율적이다.

또한, 의존 관계가 없는 규칙들의 경우 병렬로 처리될 수 있기에 수행속도가 빨라질 수 있다. 앞의 예에서 각 텍스트 파일의 단어 수를 세는 것은 서로 독립적인 작업이기에, (데이터의 크기가 크면) 다음과 같이 코어수를 늘려 동시에 실행하면 처리 속도가 빨라질 것이다.

```
snakemake --cores 3
```

### 활용 팁

앞의 예제를 좀 더 가다듬어 가면서 Snakemake 의 활용 팁을 살펴보자.

#### 기본 규칙 만들기

Snakemake 를 호출할 때마다 매번 타겟 파일을 명시하는 것은 번거롭다. Snakemake 는 타겟 파일이 없으면 첫 번째 규칙을 실행하는데, 이것을 이용해 다음처럼 기본 타겟 파일을 정의할 수 있다.

```
rule all:
    input: "temp/wc_all.png"
```

위 규칙을 `Snakefile` 의 가장 처음에 넣으면, 다음처럼 타겟 파일없이 실행할 때 자동으로 `wc_all.png` 를 타겟으로 하게 된다.

```
$ snakemake --cores 1
```

### 의존성 심화

지금까지 막연하게 의존성을 이야기했으나, 사실 Snakemake 에서 다루는 의존성은 크게 두 가지로 나누어 설명할 수 있다. 그것은 **데이터 의존성** 과 **워크플로우 의존성** 이다.

#### 데이터 의존성

데이터 의존성은 결과 파일 생성 후 입력 파일 (데이터) 이 변경되었다면 빌드가 무효화되는 것을 말한다. 이것은 결과 파일과 입력 파일의 날자를 기준으로 검사한다. Snakemake 는 데이터가 변하면 자동으로 체크하여 빌드를 다시 해준다. 예를 들어 아래와 같이 `temp/A.txt` 파일의 내용을 바꾸고,

```
Hello everybody.
```

Snakemake 를 다시 호출하면, 변경된 `data/A.txt` 에 대한 `count` 규칙만 수행되고, 그 결과인 `temp/wc_A.txt` 에 의존하는 `temp/wc_all.csv` 와 그것에 의존하는 `temp/wc_all.png` 순으로 빌드가 진행되는 것을 출력 메시지에서 볼 수 있다. 결과물인 그래프를 보아도 `A.txt` 항목의 단어수가 2 로 바뀐 것을 확인할 수 있다.

![단어수 그래프 2](/assets/2020-04-09-12-57-20.png)

#### 워크플로우 의존성

워크플로우 의존성은 `Snakefile` 에 기술된 규칙이나, 그것을 처리하는 코드에 변경이 있으면 빌드가 무효화되는 것을 말한다. Snakemake 는 워크플로우가 변하는 것은 자동으로 체크하지 못하기에, 유저가 직접 검사하여 빌드를 해주어야 한다. 다음과 같은 것이 있다:

* 규칙의 입력 변경 - `input` 의 내용이 바뀌면 무효화
  * `snakemake --list-input-changes` 로 영향 받는 출력물을 확인 가능
* 규칙의 패러미터 변경 - `params` 의 내용이 바뀌면 무효화
  * `snakemake --list-params-change` 로 영향 받는 출력물을 확인 가능
* 규칙의 스크립트 변경 - 시행 파트의 내용이 바뀌면 무효화
  * `snakemake --list-code-changes` 로 영향 받는 출력물을 확인 가능

예를 들어 `touch` 를 통해 `concat.py` 파일의 수정 날자를 갱신하고,

```
$ touch concat.py
```

Snakemake 를 다시 실행해보면, 코드가 바뀌었음에도 빌드할 것이 없다고 나온다.

```
$ snakemake --cores 1

Building DAG of jobs...
Nothing to be done.
Complete log: /Users/haje01/drills/snakemake/.snakemake/log/2020-04-13T111413.692898.snakemake.log
```

이것은 Snakemake 가 코드가 바뀐 것을 인식하지 못하기 때문이다. 다음처럼 코드 변화를 명시적으로 검사할 수 있다.

```
$ snakemake --list-code-changes

Building DAG of jobs...
temp/wc_all.png
temp/wc_all.csv
```

`concat.py` 가 변했으니, 그것의 출력인 `temp/wc_all.csv` 과 이 파일에 의존하는 `temp/wc_all.png` 가 모두 빌드되어야 한다고 알려준다.

워크플로우 내 코드 변화를 매번 명시적으로 검사하고, 필요한 경우 빌드를 해주기 위해서는 아래와 같이 명령하면 된다.

```
snakemake --cores 1 -R `snakemake --list-code-changes`
```

워크플로우 내 입력 및 매개변수 변화에 대해서도 위와 같은 방식으로 대응할 수 있다.


#### 스크립트파일 합치기

위 예제에서는 `concat` 과 `plot` 규칙을 위해 각각 파이썬 파일을 하나씩 만들어 주었다. 만약, 파이썬 파일 하나로 모든 규칙에 대응하고 싶다면 다음과 같이 하면 될 것이다. 먼저 `Snakefile` 을 다음처럼 수정하고:

```python
FILENAMES = ['A', 'B', 'C']

rule all:
    input:
        "temp/wc_all.png"

rule count:
    input:
        "data/{filename}.txt"
    output:
        "temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    input:
        expand('temp/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/wc_all.csv"
    script:
        "main.py"  # <-- 스크립트 파일명 변경

rule plot:
    input:
        "temp/wc_all.csv"
    output:
        "temp/wc_all.png"
    script:
        "main.py"  # <-- 스크립트 파일명 변경

```

`concat.py` 와 `plot.py` 의 내용을 다음처럼 `main.py` 파일로 옮긴다.

```python
import re
import pandas as pd

PTRN = re.compile(r'\s+(\d+)\sdata/(.*)$')


def concat():
    with open(snakemake.output[0], 'wt') as f:
        f.write('fname, count\n')
        for fn in snakemake.input:
            line = open(fn, 'rt').read()
            cnt, fn = PTRN.search(line).groups()
            f.write(f'{fn}, {cnt}\n')


def plot():
    df = pd.read_csv(snakemake.input[0], index_col='fname')
    plot = df.plot(kind='bar', rot=45)
    fig = plot.get_figure()
    fig.tight_layout()
    fig.savefig(snakemake.output[0])


if __name__ == '__main__':
    globals()[snakemake.rule]()  # <-- 현재 규칙에 맞는 함수 호출
```

Snakemake 를 통해 실행되는 경우 활성화되는 `snakemake` 인스턴스의 `rule` 속성에 현재 실행되는 규칙의 이름이 온다는 것을 이용하였다.

#### 특정 폴더내 모든 파일을 입력으로 하기

앞의 예에서는 어떤 입력파일 들이 있는지 알고 있다고 전제했다. `A`, `B`, `C` 가 그것이었다. 만약 확장자는 같지만 임의 파일명으로 여러 파일이 있다면 어떻게 해야할까? 파이썬의 `glob` 모듈 은 패턴에 매칭되는 파일들을 찾아주는데, Snakemake 가 제공하는 `glob_wildcards` 도 비슷한 역할을 한다.

```
FILENAMES = ['A', 'B', 'C']
```

```
FILENAMES, = glob_wildcards('temp/wc_{filename}.txt')
```

> 객체가 아닌 파일명 리스트를 얻기 위해 `FILENAMES` 다음에 `,` 가 있음에 주의하자.

#### 와일드카드에 정규식으로 제약걸기

특정 타겟을 빌드할 때 타겟 이름을 요소별로 분리해 와일드카드에 배정하고 싶을 때가 있다. 설명을 위해, 앞의 예제를 다음과 같이 바꾸어 생각해보자.

* 일별 파일들이 `YYYYMMDD` 형식 디렉토리 아래에 저장
* 일변 단어수 그래프를 `YYYYMMDD.png` 형식으로 빌드

예를 들어 다음처럼 데이터가 들어온다고 하자:

```
data/
    20200401/
        A.txt
        B.txt
        C.txt

    20200402/
        A.txt
        B.txt
        C.txt

    20200403/
        A.txt
        B.txt
        C.txt
```

2020년 4월 3일의 결과물을 `snakemake --cores 1 temp/20200403.png` 명령으로 얻으려면 `Snakefile` 의 입력은 `year`, `month`, `day` 의 세 가지  와일드카드를 이용해 다음과 같이 기술될 수 있을 것이다:

```python
FILENAMES = ['A', 'B', 'C']


rule count:
    input:
        "data/{year}{month}{day}/{filename}.txt"
    output:
        "temp/{year}{month}{day}/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    input:
        expand('temp/{{year}}{{month}}{{day}}/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/{year}{month}{day}/wc_all.csv"
    script:
        "main.py"

rule plot:
    input:
        "temp/{year}{month}{day}/wc_all.csv"
    output:
        "temp/{year}{month}{day}/wc_all.png"
    script:
        "main.py"
```

그러나, 이렇게 하면 Snakemake 는 `year` 가 4 자리, `month` 가 2 자리, `day` 가 2 자리인 것을 모르기 때문에, 잘못된 와일드카드 배정으로 인한 오류나, `RecursionError` 가 발생할 수도 있다. 따라서, 안전한 방법은 `{{와일드카드_이름, 정규식}}` 형식으로 와일드카드의 패턴을 정규식으로 지정하는 것이다. 이것을 적용하면 다음과 같다:

```python
FILENAMES = ['A', 'B', 'C']

rule count:
    input:
        "data/{year}{month}{day}/{filename}.txt"
    output:
        "temp/{year,\d{4}}{month,\d{2}}{day,\d{2}}/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    input:
        expand('temp/{{year}}{{month}}{{day}}/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/{year,\d{4}}{month,\d{2}}{day,\d{2}}/wc_all.csv"
    script:
        "main.py"

rule plot:
    input:
        "temp/{year}{month}{day}/wc_all.csv"
    output:
        "temp/{year,\d{4}}{month,\d{2}}{day,\d{2}}/wc_all.png"
    script:
        "main.py"
```

> 와일드카드의 패턴을 정규식으로 제한하는 것은 `output` 에서만 사용할 수 있다. 생각해보면 당연한 것인데, 출력 와일드카드가 정해지면 입력은 그것을 그대로 가져다 쓰면 되기 때문이다.

## TODO

* 입/출력이 하나 이상인 경우 끝에 ,
* 입/출력이 많은경우 리스트가 아닌 사전식으로
* `params` 설명
* 빌드 중 로그, 메시지 출력
* 병렬처리 추가 설명 (`threads`)
* Jupyter 노트북 실행하기
* DAG 그려보기
    * `snakemake --cores 1 wzdat-seoul/temp/ftlake/mo2/20191120.parquet --dag | dot -Tpng -o dag.png`
* S3에서 파일 입출력
  * s3에 있는 파일을 읽거나, 빌드 결과를 s3에 올릴 수 있음
  * 빌드 타겟은 `s3` 없이 로컬 파일명 형태로 지정해야 함
* 서브 워크플로우 참조
* 임시 파일 활용
* 빌드 결과물 모두 지우기
    * `$ snakemake --delete-all-output`
* 강제로 모든 규칙을 다 실행하기
* 빌드 과정 모니터링
