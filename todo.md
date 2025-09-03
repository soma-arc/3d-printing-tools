# GenMesh: GLB デフォルト出力（tinygltf）実装計画

目的: GenMesh をデフォルトで `.glb` 出力にし、`.obj` と `.vdb` はオプション指定時のみ出力。GLB 生成は syoyo/tinygltf を用いる。

## 1) CLI/仕様変更
- 既存: 常に `<basename>.vdb` と `<basename>.obj` を出力。
- 変更: 既定は `<basename>.glb` のみ出力。
- 新規フラグ:
  - `--obj` … OBJ を追加出力
  - `--vdb` … VDB を追加出力
  - （任意）`--no-glb` … GLB 出力を抑制（将来拡張、必要なら）
- `--out|-o` の basename 仕様は維持。未指定時は入力 JSON の stem を利用。
- `--help` にフラグ説明を追記。

## 2) 依存とビルド設定（Premake）
- 依存: `GenMesh/src/tiny_gltf.h`（追加済）。リンクは不要（ヘッダオンリー）。
- 画像は扱わないため、STB は不要だが、`stb_image*.h` も同梱済のためどちらでも可。
  - STB を使わない場合: `defines { "TINYGLTF_NO_STB_IMAGE", "TINYGLTF_NO_STB_IMAGE_WRITE" }` を追加。
  - STB を使う場合: 定義は不要（今回の GLB 書き出しのみならどちらでもビルド可能）。

## 3) GLB エクスポータの追加（tinygltf）
- 新規関数: `writeGlb(const std::string& basename, const std::vector<openvdb::Vec3s>& points, const std::vector<openvdb::Vec3s>& normals, const std::vector<openvdb::Vec4I>& quads, const std::vector<openvdb::Vec3I>& triangles)`。
- 実装方針:
  1. 面の三角化: `quads` → `(x,y,z)` と `(x,z,w)` の 2 三角面へ展開。`triangles` はそのまま結合。
  2. バッファ生成:
     - positions: `float32 * 3 * points.size()`
     - normals:   `float32 * 3 * normals.size()`（未生成時は省略可だが現状は生成あり）
     - indices:   `uint32 * (三角面数*3)`
  3. tinygltf の `Model/Buffer/BufferView/Accessor/Mesh/Primitive/Node/Scene` を構築。
     - Accessor 設定:
       - POSITION: `type: VEC3`, `componentType: FLOAT`, `count: points.size()`, `min/max` を AABB から算出。
       - NORMAL:   `type: VEC3`, `componentType: FLOAT`, `count: normals.size()`
       - INDICES:  `type: SCALAR`, `componentType: UNSIGNED_INT`, `count: indexCount`
     - BufferView `target`: vertices → `ARRAY_BUFFER (34962)`, indices → `ELEMENT_ARRAY_BUFFER (34963)`。
     - Primitive: `attributes["POSITION"]`, `attributes["NORMAL"]`, `indices`、`mode = TRIANGLES (4)`。
  4. 書き出し: `TinyGLTF gltf;` → `gltf.WriteGltfSceneToFile(model, basename + ".glb", /*embedImages=*/true, /*embedBuffers=*/true, /*pretty=*/false, /*writeBinary=*/true)` を使用。

## 4) 既存フローへの統合
- `main.cpp` の引数定義に `args::Flag exportObj`, `args::Flag exportVdb` を追加。
- `makeMesh(...)` 終了部の出力を条件分岐:
  - 既定: `writeGlb(outputBasename, points, normals, quads, triangles);`
  - `if (exportObj) writeObj(...);`
  - `if (exportVdb) writeGrid(...);`
- `--help` のメッセージ更新。

## 5) ドキュメント更新
- `GenMesh/README.md`: 既定が GLB 出力である旨、利用例を追記。
- 追加出力のフラグ例（`--obj`, `--vdb`）を記載。

## 6) 動作確認
- ビルド: `cd GenMesh && premake5 gmake && make -C build config=Release`。
- 実行: 既存 JSON を入力に
  - 既定 → `<basename>.glb` のみ生成
  - `--obj` → `.glb` と `.obj`
  - `--vdb` → `.glb` と `.vdb`
- GLB 検証: ファイル種別/サイズ確認、外部の glTF Viewer/Validator で読み込み確認。

## 7) 実装注意
- インデックスは 0 起点（glTF 準拠）。
- バッファの 4 バイト境界に注意（`byteOffset` は 4 の倍数）。
- 法線が空の場合は `NORMAL` 自体を省略可。現状は `computeNormals` 出力を使用。
- 例外/失敗時は標準エラー出力し、終了コードは既存ポリシーに合わせる。
