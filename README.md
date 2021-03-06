# K-dimensional tree space-partitioning data structure

## Overview

The k-d tree is a binary tree in which every node is a k-dimensional point. Every non-leaf node can be thought of as implicitly generating a splitting hyperplane that divides the space into two parts, known as half-spaces. Points to the left of this hyperplane are represented by the left subtree of that node and points right of the hyperplane are represented by the right subtree.

Read more on [wikipedia](https://en.wikipedia.org/wiki/K-d_tree)

Additional specific features are implemented in the current library for collision detection tasks. By splitting objects bounding volumes it allows not only to detect the fact of collision, but also to locate a collision point for complex objects with the given accuracy.

Visualizer: [kdtree-demo](https://github.com/swizard0/kdtree-demo)

![kdtree visualizer](https://github.com/swizard0/kdtree-demo/blob/master/images/screenshot_00.png "KD-tree visualizer")

## Example usage

See library [tests](src/lib.rs#L365) section or [kdtree demo](https://github.com/swizard0/kdtree-demo) visualizer code.

## Authors

* **Alexey Voznyuk** - *Initial work* - [swizard0](https://github.com/swizard0)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
