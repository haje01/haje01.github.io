---
layout: post
title: "Argo 기반 ETL 코드에 Skaffold 적용하기"
description: 
date: 2023-03-27
tags: [draft]
---
{% raw %}
이 글은 앞서 소개한 [Argo Events 와 Workflows 를 이용한 ETL]() 의 예제에 Skaffold를 적용한 것에 대해 설명한다. 따라서 그 문서를 꼭 먼저 읽고 이 길을 읽기 바란다.

## 준비 

이 글은 Linux 환경을 전제로 설명한다.

## Skaffold 설치 

[Skaffold 공식 페이지](https://skaffold.dev/) 를 참고하여 설치하자.

### Skaffold 를 사용해야 하는 이유

[Skaffold](https://skaffold.dev/) 는 컨테이너 및 쿠버네티스 기반의 지속적 개발을 제공하는 CLI 툴이다. 크게 다음과 같은 종류의 일을 해준다:
- 빌드
  - 필요한 경우 컨테이너 이미지를 빌드하고 푸쉬해준다.
  - 파일 복사로 충분한 경우 실행중인 컨테이너에 복사해준다.
- 배포 - 쿠버네티스 리소스를 설치한다.
  - `kubectl` 을 통해 YAML 매티페스트 파일을 설치
  - `helm`, `kustomize` 를 통한 설치
- 후킹
  - 빌드 및 배포의 시점에 특정 작업을 끼워넣을 수 있다.
  - 스크립트 등을 통한 전/후 처리에 유용하다.

Skaffold 는 대상 파일을 모니터링하여 파일 변경시 적절한 빌드나 배포를 자동으로 진행 (소위 Hot-reload 기능) 해 개발시 많은 도움이 된다:

- 상황별 초기 설치작업을 자동화
- Dockerfile 에 변경이 있을 때 자동으로 빌드해주고, 그것을 사용하는 리소스를 재생성 
- 소스 코드가 변경되면 자동으로 그것을 이용하는 컨테이너에 복사

Skaffold 의 이런 기능을 이용해서 개발을 하면서 실시간으로 리소스가 자동으로 갱신되는 것을 **지속적 개발 (Continuous Developing)** 모드라 하는데, 처음부터 쿠버네티스 환경에서 개발을 시작하는 경우 없어서는 안될만큼 편리하다.

혹자는 *꼭 쿠버네티스 환경 안에서 개발을 할 필요는 없지 않나?* 라고 생각할 수도 있겠다. 일반적인 환경에서 개발하고 배포를 위해 쿠버네티스에 맞게 패키징하는 방식을 선택할 수도 있다. 그렇지만, 프로젝트 초기부터 쿠버네티스 배포를 상정하는 경우가 늘고 있고, 라이브 후에도 지속적인 기능 추가와 디버깅까지 고려하면 쿠버네티스 환경 안에서의 개발에 익숙해지는 것도 좋을 것이다.

> 컨테이너 환경에서 최신 소스 코드를 확보하기 위해서는 여러가지 방법이 있다. 대표적으로는 다음과 같다:
>
> - 컨테이너 이미지에 포함
> - 호스트의 경로를 마운팅
> - 컨테이너 시작시 최신 코드를 git 으로 받기
> 
> 호스트의 경로를 마운팅하는 식은 해당 노드에서만 코드에 접근할 수 있기에 멀티노드에서는 사용할 수 없다. 컨테이너 시작시 최신 코드를 받는 방법은 제한적으로 이용할 수도 있으나 다음과 같은 문제점이 있다:
> - 컨테이너의 장점인 빠른 시작 시간이 느려지는 점
> - 같은 이미지 버전을 사용하는 컨테이너들도 배포 시점에 따라 서로 다른 코드가 동작할 가능성
> - git 관련 로그인 정보가 컨테이너에서 노출될 가능성 
> 
> 따라서, 컨테이너 이미지를 만들때 코드를 함게 포함하여 빌드하고 태그를 붙여주는 방식이 가장 추천된다. 

이 글은 이전 문서의 예제에 Skaffold 를 적용한 것을 주로 설명한다. Skaffold 에 관한 모든 내용을 설명하는 것은 이 글의 범위를 넘어서는 것이기에, 자세한 것은 공식 페이지와 [관련 블로그](https://blog.slashuniverse.com/25) 를 참고하도록 하자. 

## 예제 설치 

```sh
git clone https://github.com/haje01/argo-etl.git
```

예제의 폴더 구조는 아래와 같다.
```
argo-etl/
  # 인프라 셋업 
  setup/  
  # 개별 예제들
  minio/
  minio-parq/
  s3/
  ftp/
```

`argo-etl` 아래 최상위 폴더들은 각각의 Skaffold 프로젝트를 위한 것이다. 단, `setup` 프로젝트의 경우 다른 예제 프로젝트에서 필요한 초기화를 대신해준다.

### 프로젝트 폴더 구조 

`setup` 프로젝트 폴더를 예로 Skaffold 프로젝트의 구조를 한 번 살펴보자.

```
setup/
  skaffold.yaml
  k8s/ 
  local/
  eks/
```

`skaffold.yaml` 파일은 Skaffold 의 프로젝트 설정 파일이다. 그 안에는 빌드, 테스트, 앱 배포에 관련된 내용을 기술할 수 있다.

쿠버네티스 관련 매니페스트 파일은 관례적으로 `k8s/` 폴더에 둔다. 컨테이너 별로 필요한 파일은 `local` 이나 `eks` 처럼 별도 폴더를 만들어서 저장하는데, 이것을 **컨피그 (Config) 폴더** 라 한다. 컨피그 폴더는 필요에 따라 생성해 이용하면 된다.

`setup` 프로젝트의 경우 클러스터 생성 후 설정에 관한 프로젝트이다. 로컬 클러스터 설정에 필요한 파일은 `local` 컨피그 폴더에, Amazon EKS 클러스터 설정에 필요한 파일은 `eks` 컨피그 폴더에 있으며, 공통적으로 필요한 매니페스트 파일은 `k8s` 폴더에 있다.

### 프로젝트 파일의 구조 

`skaffold.yaml` 파일을 보통 다음과 같은 구조를 가진다.

```yaml
apiVersion: skaffold/v4beta3
kind: Config
build:
  # 생성할 이미지 등 빌드 관련 내용 
manifests:
  # 대상 매니페스트 파일 
deploy:
  # 배포 관련 내용 
```

Skaffold 는 대상 매니페스트 파일을 쿠버네티스 클러스터로 보내기 전에 **렌더링 (Rendering)** 을 한다. `manifest` 섹션에는 대상 매니페스트 파일을 기술하고, 그것을 어떻게 렌더링할지는 `deploy` 에서 기술한다.

상황별로 빌드, 테스트, 배포 설정을 다르게 기술하고 싶은 경우는 **프로파일 (Profile)** 을 이용한다. `setup` 프로젝트의 `skaffold.yaml` 을 보면 다음처럼 `local` 및 `eks` 프로파일로 나누어져 있다.

```yaml
apiVersion: skaffold/v4beta3
kind: Config
profiles:
  # 로컬 (Minikube) 인프라 설치
- name: local
  manifests:
  deploy:
  # AWS EKS 인프라 설치
- name: eks
  manifests:
  deploy:
```

`local` 에는 로컬 클러스터용 설정, `eks` 에는  AWS EKS 클러스터용 설정이 있는데, 프로파일별로 `manifests` 와 `deploy` 를 다르게 지정하고 있는 것을 알 수 있다.

> `setup` 프로젝트의 경우 프로파일과 컨피그 폴더의 형태가 일치하나, 항상 이런 것은 아니다.

컨피그 폴더에는 특정 프로파일에만 필요한 쿠버네티스 매니페스트 파일도 위치한다. 예를 들어 이전 문서에서 설명한 MinIO 를 위한 아티팩트 리소스 설정은 `setup/local/art-repo.yaml` 에 있다.

## 로컬 클러스터

minikube 를 이용해 로컬 클러스터가 생성되어 있다는 가정아래, 클러스터를 설정하고 관련 예제들을 살펴보자.

### 로컬 클러스터 설정

먼저 필요한 패키지 설치 등 설정을 진행해야한다. `setup` 프로젝트의 `local` 컨피그 폴더에 로컬 클러스터 환경을 위한 공통 설정이 기술되어 있다. `skaffold.yaml` 에는 로컬 `local` 및 AWS EKS 용 `eks` 프로파일을 구분하여 기술하고 있는데, `-p local` 식으로 특정 프로파일을 선택할 수 있다.

`argo-etl/setup` 폴더로 이동 후 다음처럼 실행하면 된다.

```
skaffold dev -p local
```

> 만약 실행중 에러가 발생하면 Skaffold 는 즉시 실행을 정지하고, 그때까지 설치된 내용을 삭제한다.

다음과 같은 작업이 진행된다:

- Argo Events 설치 
- Argo Events 의 커스텀 리소스인 EventBus 설치 
- Argo Events 를 위한 서비스 어카운트 (`operate-workflow-sa`) 생성
- Argo Workflows 설치 
- Argo Workflows UI 토큰을 얻기 위한 Secret 생성
- MinIO 설치
- Argo Workflows 의 기본 아티팩트 저장소 (여기서는 MinIO) 설정

Skaffold `dev` 명령에서는 Control + C 키로 명령을 종료하면 설치된 모든 내용도 함께 삭제가 되는 점에 주의하자.

### MinIO 예제 

`minio` 폴더에는 MinIO 의 특정 버킷 경로에 입력 CSV 파일이 올라오면 ETL 후 CSV 파일로 결과를 저장하는 프로젝트가 있다. 

`skaffold.yaml` 파일의 내용은 다음과 같다. 

```yaml
```

`minio` 폴더로 이동 후 아래와 같이 Skaffold 를 호출하면,

```
skaffold dev
```

이전 문서에서 수동으로 했던 다음과 같은 작업이 자동으로 진행된다:

- 대상 MinIO 버킷에 .CSV 파일이 올라 오는 것을 검출하기 위한 이벤트 소스 `minio` 생성 
- 이벤트가 발생하면 ETL 워크플로우를 호출하기 위한 센서 `minio` 생성
- `etl.py` 를 실행하기 위한 컨테이너 이미지 `minio` 빌드
- MinIO 에 예제용 버킷 `etlproj` 생성 

정상적으로 작업이 완료되면, 이전 문서와 같은 식으로 MinIO UI 를 통해 CSV 파일을 업로드해보자. ETL 결과가 CSV 파일로 남는 것을 확인할 수 있을 것이다.

### MinIO Parquet 예제

`minio-parq` 폴더에는 MinIO 의 특정 버킷 경로에 입력 CSV 파일이 올라오면 ETL 후 Parquet 파일로 결과를 저장하는 프로젝트가 있다. 

`skaffold.yaml` 파일의 내용은 다음과 같다. 

```yaml
```

폴더로 이동 후 아래와 같이 Skaffold 를 호출하면 예제를 위한 초기화 작업이 진행된다. 

```
skaffold dev
```

이전 문서에서 수동으로 했던 다음과 같은 작업이 진행된다:

- 대상 MinIO 버킷에 .CSV 파일이 올라 오는 것을 검출하기 위한 이벤트 소스 `minio-parq` 생성 
- 이벤트가 발생하면 ETL 워크플로우를 호출하기 위한 센서 `minio-parq` 생성
- `etl.py` 를 실행하기 위한 컨테이너 이미지 `minio-parq` 빌드
- MinIO 에 예제용 버킷 `etlproj` 생성 

정상적으로 작업이 완료되면, 이전 문서와 같은 식으로 CSV 파일을 업로드해 ETL 결과가 Parquet 파일로 남는 것을 확인할 수 있을 것이다.

## Amazon EKS 클러스터

단일 EC2 인스턴스에 쿠버네티스 클러스터를 만들어 진행했던 이전 문서와 달리, 이번에는 AWS 제공하는 관리형 쿠버네티스 서비스인 [Amazon Elastic Kubernetes Service](https://aws.amazon.com/ko/eks/) (줄여서 EKS) 클러스터를 이용하여 `s3` 및 `ftp` 예제를 살펴보겠다. 

이를 위해 먼저 [AWS CLI](https://docs.aws.amazon.com/ko_kr/cli/latest/userguide/getting-started-install.html) 및 [eksctl](https://eksctl.io/introduction/#installation) 의 설치가 필요하다. 

링크를 참고하여 설치하고 준비가 되었으면 EKS 클러스터 생성부터 시작하겠다.
 
### EKS 클러스터 생성

다음처럼 클러스터를 생성하자.

```sh
eksctl create cluster --name=prod --nodes=1 --node-type m5.large --node-volume-size=20
```

`m5.large` 타입에 디스크 20GB 를 가지는 EC2 인스턴스 (노드) 하나로 구성된 `prod` 이라는 이름의 클러스터를 생성하라는 명령이다. EKS 클러스터는 초기 생성에 시간이 꽤 걸리는데, 15 분 정도 기다리면 클러스터가 만들어질 것이다.

> EKS 클러스터의 이용이 끝나면 비용을 절감하기 위해 꼭 삭제해주도록 하자. 일반적으로는 `eksctl` 을 사용해서 다음과 같이 삭제한다.
> ```sh
> eksctl delete cluster prod
> ```
> 그런데, `eksctl` 은 클러스터 삭제에 실패하거나, 삭제에 성공하더라도 VPC 관련 리소스가 남아 있는 경우가 있다. 이런 경우 `setup/eks-delete.sh` 를 통해 정리할 수 있다. 다음과 같이 호출하면 된다.
> ```
> EKS_CLUSTER=prod sh eks-clean.sh
> ```
> 위 스크립트도 완벽하지 않아 원하지 않는 결과가 나올 수도 있기에, 중요한 리소스에 대해서는 하나씩 찾아 수작업으로 지워야 한다.

제대로 설치가 되었으면 아래 명령으로 EKS 클러스터가 추가되었음을 알 수 있고,

```sh
$ kubectl config get-clusters

NAME
minikube
prod.ap-northeast-2.eksctl.io
```

기존 minikube 대신 EKS 클러스터가 활성화 되었음을 확인할 수 있다.

```sh
$ kubectl config current-context 

<IAM 사용자 이름>@prod.ap-northeast-2.eksctl.io
```

AWS UI 의 EKS 페이지에서도 클러스터를 확인할 수 있다. 

![EKS Cluster](/assets/2023-03-29-11-43-29.png)

다음은 자신의 AWS 계정 ID (숫자로 된 것) 를 환경변수로 노출하고,

```sh
export AWS_ACCOUNT_ID=<AWS 계정 ID>
```

Skaffold 의 기본 컨테이너 이미지 저장소를 AWS ECR 로 지정하자. 

```sh
skaffold config set default-repo $AWS_ACCOUNT_ID.dkr.ecr.ap-northeast-2.amazonaws.com
```

Skaffold 를 통한 컨테이너 이미지 빌드시 ECR 에 Push 하기 위해서는 아래와 같은 로그인도 필요하다.

```
aws ecr get-login-password | docker login --username AWS --password-stdin $AWS_ACCOUNT_ID.dkr.ecr.ap-northeast-2.amazonaws.com
```

### EKS 클러스터 셋업

추가적으로 EKS 클러스터에 설치 및 설정을 진행해야 하는데, 편리하게 `setup` 프로젝트의 `eks` 컨피그 폴더에 EKS 클러스터 환경을 위한 공통 설정이 기술되어 있다. `argo-etl/setup` 폴더로 이동 후 다음처럼 명령하면,

```
skaffold dev -p eks
```

기술된 내용이 실행되게 된다. 

위 과정에서는 크게 두 가지 종류의 설정이 진행되는데 하나는 EKS 클러스터 관련 설정이고, 다른 하나는 패키지 설치이다.

EKS 클러스터 설정은 `eks/init.sh` 파일에 기술되어 있는데, 다음과 같은일이 진행된다:

1. OIDC 프로바이더 연동
2. `AWSLoadBalancerControllerIAMPolicy` 생성
3. `AmazonEKSLoadBalancerControllerRole` 역할을 만들고, 그것과 연결된 쿠버네티스 서비스 어카운트 생성

1 번은 클러스터에 대해 한 번만 수행되고, 2 번과 3번은 AWS 계정에 대해 한 번만 수행되면 된다. 

> 기존에 작업된 내용이 있으면 현재 작업은 스킵된다. 

이 부분이 로컬에서 minikube 를 사용하는 것에 비해 복잡한 과정이나, Ingress 및 Scale In/Out 지원 등 실제 서비스를 운영할 수 있을 정도로 강력한 쿠버네티스 환경임을 감안하면 납득이 될 것이다.

다음으로 진행되는 `eks` 폴더 아래의 Helm 패키지 설치는 로컬 클러스터의 경우와 크게 다르지 않다. 다음과 같은 작업이 진행된다:

- AWS 로드밸런서 컨트롤러 (`aws-load-balancer-controller`) 설치 
- Argo Events 설치 
- Argo Events 의 EventBus 설치 
- Argo Events 를 위한 ServiceAccount (`operate-workflow-sa`) 생성
- Argo Workflows 설치 
- Argo Workflows UI 토큰을 얻기 위한 Secret 생성

Argo Workflows UI 페이지는 아래와 같이 인그레스의 주소를 얻어서,

```
$ kubectl get ingress
NAME                        CLASS    HOSTS                                ADDRESS                                                             PORTS   AGE
awf-argo-workflows-server   <none>   *.ap-northeast-2.elb.amazonaws.com   k8s-public-a25c3d33fb-1538911396.ap-northeast-2.elb.amazonaws.com   80      2m31s
```

다음과 같이 토큰을 얻어 로그인하면 된다.
```
SECRET=$(kubectl get sa awf-argo-workflows-server -o=jsonpath='{.secrets[0].name}')
ARGO_TOKEN="Bearer $(kubectl get secret $SECRET -o=jsonpath='{.data.token}' | base64 -d)"
echo "$ARGO_TOKEN"
```

Control + C 키로 `dev` 명령을 종료하면 설치된 패키지가 삭제된다.

> 만약 라이브 배포 목적 등의 이유로 설정 패키지가 유지되기를 원한다면, 다음처럼 `deploy` 를 이용하자.
> ```
> skaffold deploy -p eks
> ```
> `deploy` 로 생성된 리소스는 `delete` 로 삭제할 수 있다.
> ```
> skaffold delete -p eks
> ```

### S3 예제

`s3` 폴더에는 S3 특정 버킷 경로에 CSV 파일이 올라오면 ETL 후 S3 에 Parquet 파일로 저장하는 프로젝트가 있다. 이전 문서의 예제와 거의 같으나 EKS 클러스터와 Skaffold 를 이용한다는 점이 다르다.

ETL 위한 파이썬 컨테이너 이미지를 빌드 후 푸쉬하게 되는데, 이를 위해 ECR 에 `s3-etl` 이라는 저장소 이름 (Repository) 을 프라이빗으로 만들어 둔다.

> 저장소 관련 용어:
> Registry (레지스트리) - 저장소
> Repostiroy (리포지토리) - 저장소 이름

![ECR Repository](/assets/2023-04-04-16-57-02.png)

예제에서 파일을 올리고 받기 위해 S3 버킷이 필요하다. 중복되지 않는 적당한 이름으로 미리 만든다. 

`skaffold.yaml` 파일의 내용은 아래와 같다.

```yaml
apiVersion: skaffold/v4beta3
kind: Config
metadata:
  name: s3
build:
  # ETL 코드를 실행할 컨테이너 이미지. `etl` 폴더에 정보가 있다.
  artifacts:
  - image: s3-etl
    context: etl
deploy:
  # 단순 YAML 파일 복사가 아닌, Helm 차트를 통해 배포한다.
  helm:
    # 패키지 설치전 처리
    hooks: 
      before:
      - host:
          command: ["sh", "init.sh"]
          os: [darwin, linux]  
    releases:
    - name: s3-sns 
      chartPath: setup
      # 무한 배포를 막기 위해 변수별 파일 이용 
      valuesFiles:
      - vals/topic-arn.yaml
      - vals/ingress-addr.yaml
      # 빌드 과정에서 생성된 컨테이너 이미지 이름을 변수로 전달한다.
      setValueTemplates:
        image: '{{.IMAGE_FULLY_QUALIFIED_s3_etl}}'
```

좀 복잡한데, 다음과 같은 일을 한다:
- ETL 코드를 실행할 컨테이너 이미지 `s3-etl` 을 빌드
- 설치 전 초기화 스크립트 `init.sh` 실행
- Helm 차트를 이용해 다음과 같은 리소스를 설치
  - Argo Events 의 이벤트 소스 
  - Argo Events 의 센서 (실행할 워크플로우 포함)
  - 인그레스 

`init.sh` 스크립트는 이전 문서에서 소개한 복잡한 AWS S3 및 SNS 설정 과정을 자동화해준다. 

배포를 단순 YAML 형 매니페스트 파일이 아닌 Helm 차트를 이용하는 것은, 초기화 스크립트에서 생성한 리소스를 매니페스트에서 활용하기 위함이다.

초기화 스크립트를 위해 Skaffold 호출 전 다음과 같은 환경변수를 결정해야 한다:

- `PROJECT` - 구분할 수 있는 프로젝트 명 
- `S3_BUCKET` - 대상 (입력/출력) S3 버킷 명
- `S3_PREFIX` - 대상 S3 입력 객체 접두사

설명에서는 프로젝트 명은 `myproj` 로 하고, `my-bucket` 버킷의 `input/` 폴더 아래에 `.csv` 파일이 올라오는 것으로 가정하겠다. `s3` 폴더에서 아래와 같이 환경 변수를 지정하면서 Skaffold 를 호출하면 된다. 

```
PROJECT=myproj S3_BUCKET=my-bucket S3_PREFIX=input/ skaffold dev
```

> 만약 ECR 로그인이 되어있지 않거나, 저장소 이름이 없으면 이미지 Push 과정에서 에러가 발생할 수 있다. 앞의 내용을 참고하여 진행하도록 하자.

`skaffold dev` 를 호출하면 `init.sh` 스크립트가 다음과 같은 일을 해준다:

1. SNS 토픽 `myproj-s3-noti` 생성
2. 토픽 정책 `myproj-s3-noti-policy` 설정 
3. 대상 버킷에 알림 구성
4. 알림 전용 IAM 유저 `myproj-noti-bot` 생성
5. 유저가 SNS 토픽 알림을 받을 수 있도록 정책 `myproj-noti-bot-policy` 생성 및 적용
6. 추가적으로 Helm 차트가 설치될 수 있도록 토픽 ARN 및 IAM 유저 키정보를 `setup/secret-values.yaml` 파일에 저장

구체적인 내용은 `init.sh` 파일을 참고하자. 

또한, Helm 차트를 통해 다음과 같은 리소스가 생성된다:

- 유저가 만든 예제용 S3 버킷에 CSV 파일 생성시 SNS 웹훅 호출을 받기위한 EventSource `s3` 생성
- 이벤트가 발생하면 ETL 워크플로우를 호출하기 위한 Sensor `s3` 생성
- `etl.py` 를 실행하기 위한 컨테이너 이미지 `s3` 빌드

정상적으로 작업이 완료되면, 이전 문서와 같은 식으로 CSV 파일을 업로드해 ETL 결과가 CSV 로 파일로 남는 것을 확인할 수 있을 것이다.

Control + C 키로 `dev` 명령을 종료하면 Helm 으로 설치된 패키지가 삭제된다.

그러나 `init.sh` 스크립트를 통해 생성된 AWS 리소스 및 설정은 제거되지 않기에, 다음과 같은 환경 변수를 지정해 삭제 스크립트를 불러주어야 한다.

```
PROJECT=myproj S3_BUCKET=my-bucket sh delete.sh
```

> 아직 Skaffold 는 [종료 (Cleanup) 훅을 지원하지 않기에](https://github.com/GoogleContainerTools/skaffold/issues/8308) 이 부분을 자동화하는 것이 어렵다.

### FTP 예제

`ftp` 폴더에는 FTP 서버 내 특정 디렉토리에 파일이 올라오면 ETL 후 S3 에 Parquet 파일로 저장하는 프로젝트가 있다. 이전 문서의 예제와 거의 같으나 EKS 클러스터와 Skaffold 를 이용한다는 점이 다르다.

ETL 위한 파이썬 컨테이너 이미지를 빌드 후 푸쉬하게 되는데, 이를 위해 ECR 에 `ftp-etl` 이라는 저장소 이름 (Repository) 을 만들어 둔다.

이 예제에서도 S3 버킷이 필요한데, 앞에서 만들어 둔것을 이용하겠다.

이전 문서에서 소개한 것처럼 FTP 서버 및 커스텀 이벤트 서버의 설치가 필요한데, 이 것도 Skaffold 프로젝트에 의해 자동화 되어있다. 다음과 같은 환경변수를 건내주어야 한다:

- `PROJECT` - 구분할 수 있는 프로젝트 명 
- `S3_BUCKET` - 결과를 저장할 S3 버킷 명
- `S3_PREFIX` - 결과를 저장할 S3 객체 접두사

`ftp` 폴더에서 아래와 같이 환경 변수를 지정하면서 Skaffold 를 호출하면 된다. 

```
PROJECT=myproj S3_BUCKET=my-bucket S3_PREFIX=output/ skaffold dev
```

다음과 같은 작업이 진행된다:

- FTP 서버 (vsftp) 설치 
- FTP 용 커스텀 이벤트서버 설치 
- Argo Events 이벤트소스 설치 
- Argo Workflows 워크플로우를 포함한 Argo Events 센서 설치 

FTP 서버 설치가 완료되면 다음처럼 접속 명령이 출력된다. 

```
## FTP 접속 명령 ##
export FTP_IP=$(kubectl get svc --namespace {{ .Release.Namespace }} -l "app.kubernetes.io/name=vsftpd" -o jsonpath="{.items[0].status.loadBalancer.ingress[0].hostname}")
ftp -p $FTP_IP
```

아래 내용으로 `test.csv` 파일을 만든 뒤,

```
1,aaa,100
2,bbb,80
3,ccc,90
```

FTP 접속 명령을 복사해 실행하면 유저 이름 `admin`, 암호 `djemals` 로 접속이 가능하다. 

> `Name or service not known` 에러가 발생하거나 동작이 멈춰있으면, 잠시 후 다시 시도해보자.


FTP 클라이언트에서 다음처럼 `test.csv` 파일을 FTP 에 올리면,

```
put test.csv
```

다음과 같은 순서로 진행이 된다.

1. 이벤트 서버가 FTP 에서 변화를 감지하면, 이벤트 소스에 전달 
2. 이벤트 소스를 통해 센서 의존성이 해결
3. 센서가 워크플로우를 생성
4. 워크플로우에서 ETL 진행

문제없이 진행되면 예제의 경우 S3 `my-bucket` 아래에서- `output/test.parquet` 파일을 확인할 수 있을 것이다. 

## 기타

### Skaffold 를 라이브 Setup 에 이용하기

지금까지는 `setup` 프로젝트에서 `skaffold dev`를 호출해 공통 설정을 진행하였다. `dev` 는 지속적 개발을 위한 명령으로 종료되면 자신이 초기화한 리소스를 제거하는 특징이 있다. 

Skaffold 에는 `deploy` 와 `delete` 라는 명령이 있는데, 이것을 이용하면 라이브 Setup 에도 활용이 가능하다.

그렇지만 사이드이펙트가 있기에 주의해서 사용해야 할 것

{% endraw %}